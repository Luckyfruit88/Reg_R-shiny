# Start from this directory: Rscript -e 'shiny::runApp(".", host="127.0.0.1", port=3838)'
needed <- c("shiny", "bslib", "DT", "data.table", "processx", "future", "promises", "jsonlite")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("缺少 R 包：", paste(missing, collapse = ", "), "。请按 README 安装。")
if (utils::packageVersion("shiny") < "1.8.1") stop("需要 shiny >= 1.8.1。")
library(shiny)
library(bslib)
backend_file <- normalizePath("backend.R", mustWork = TRUE)
backend <- new.env(parent = globalenv())
sys.source(backend_file, envir = backend)
# One app process / trusted local users. Provision a separate worker pool for production.
future::plan(future::multisession, workers = 2)
onStop(function() future::plan(future::sequential))

bam_root <- Sys.getenv("REGTOOLS_BAM_DIR", unset = file.path(getwd(), "data"))
if (!dir.exists(bam_root)) stop("REGTOOLS_BAM_DIR 不存在：", bam_root)
bam_root <- normalizePath(bam_root, mustWork = TRUE)
bams <- list.files(bam_root, pattern = "\\.bam$", ignore.case = TRUE, full.names = TRUE)
if (length(bams)) {
  bams <- normalizePath(bams, mustWork = TRUE)
  bams <- bams[startsWith(bams, paste0(bam_root, .Platform$file.sep))]
}
bam_choices <- stats::setNames(bams, basename(bams))

ui <- page_sidebar(
  title = "RegTools · Junction Explorer",
  sidebar = sidebar(width = 345,
    h5("1 · 数据与固定统计区间"),
    selectInput("source", "BAM 来源", c("内置模拟数据" = "demo", "服务器目录" = "server")),
    conditionalPanel("input.source == 'server'",
      selectInput("bam", "选择 BAM（只读）", choices = bam_choices),
      helpText("索引放在 BAM 旁边；目录由 REGTOOLS_BAM_DIR 指定。")),
    textInput("chrom", "染色体", "chrDemo"),
    numericInput("roi_start", "统计区间起点（1-based）", 101, min = 1),
    numericInput("roi_end", "统计区间终点（含该碱基）", 450, min = 1),
    tags$details(tags$summary("过滤与 RegTools 参数"),
      numericInput("mapq", "最低 MAPQ", 20, min = 0, max = 255),
      numericInput("baseq", "最低 BaseQ（只影响 depth）", 0, min = 0, max = 93),
      numericInput("anchor", "RegTools minimum anchor (bp)", 8, min = 1),
      numericInput("min_intron", "最短内含子 (bp)", 70, min = 1),
      numericInput("max_intron", "最长内含子 (bp)", 500000, min = 1),
      selectInput("strand_mode", "RegTools 文库方向", c("XS 标签" = "XS", "RF / first-strand" = "RF", "FR / second-strand" = "FR")),
      checkboxInput("exclude_dup", "排除已标记的 duplicate", FALSE),
      checkboxInput("nh1", "只保留 NH==1（缺失 NH 也排除）", FALSE),
      helpText("固定排除 unmapped、secondary、supplementary、QC-fail。MAPQ 阈值不是唯一比对的通用替代。")),
    input_task_button("run", "读取区间 / 应用过滤参数"),
    hr(),
    h5("2 · 在已读取结果上调整"),
    numericInput("target_start", "目标内含子第一个碱基（1-based）", 201, min = 1),
    numericInput("target_end", "目标内含子最后一个碱基（1-based）", 300, min = 1),
    sliderInput("delta", "Near：两端分别允许偏移 ±bp", min = 0, max = 50, value = 5, step = 1),
    helpText("Near 排除 exact；不是“任意一个剪接位点附近”。点击 junction 表的一行可设为目标。"),
    uiOutput("display_control")
  ),
  uiOutput("status"),
  uiOutput("metrics"),
  p(textOutput("denominator_note")),
  navset_card_tab(
    nav_panel("图形",
      card(card_header("Read depth · 全部位置（包括零覆盖）"), plotOutput("depth_plot", height = 250)),
      card(card_header("Junction arcs · 当前视窗内支持数最高的 30 条"), plotOutput("junction_plot", height = 260)),
      p("图形缩放只改变显示；固定统计区间与归一化分母不变。弧线不是转录本注释。")),
    nav_panel("Junction 表", DT::DTOutput("junction_table")),
    nav_panel("Read 证据",
      checkboxInput("evidence_all", "显示所有通过 RegTools 的 junction 证据（否则只显示 exact / near）", FALSE),
      DT::DTOutput("evidence_table"),
      p("每行是一条 alignment 的一个 N 操作；同一 read_id 可能多行。导出不包含 SEQ/QUAL。")),
    nav_panel("导出与核验",
      downloadButton("download_summary", "指标 CSV"),
      downloadButton("download_junctions", "Junction CSV"),
      downloadButton("download_depth", "Depth TSV"),
      downloadButton("download_evidence", "全部 read 证据 CSV"),
      downloadButton("download_metadata", "参数 / 日志 JSON"),
      hr(), verbatimTextOutput("audit"),
      h5("统计口径"),
      p("Reads 按保留的 primary alignment 记录计数，paired-end 两端分别计算；不是 fragment 或 UMI 数。所有链（包括未知链）合并汇总，junction 表保留链信息。"),
      p("分母仅计 M、=、X 比对块真正覆盖固定统计区间的 reads；只以 N 或 D 跨越区间不算覆盖。"),
      p("junction_reads 对通过 RegTools 的 junction 取 read_id 并集。exact / near 分别取匹配事件的 read_id 并集；near 不包括 exact 事件。"),
      p("每 100 reads = 100 × exact count / denominator_reads；分母为 0 时是 NA。它不是 PSI。"),
      p("Mean depth 是固定区间逐碱基 depth 的平均值（含内含子和零覆盖），另导出目标两侧最多 50 bp 的均值。")),
    nav_panel("使用边界",
      p("这是研究用局部分析原型，不是经过临床验证的软件。默认不上传 BAM、不修改原文件。"),
      p("只支持本机或可信内网使用；未实现登录、跨用户权限、任务取消和共享任务队列。不要直接暴露到公网。"),
      p("单次最多 250 kb、100,000 条候选 alignment；超过上限报错而不是抽样。目标变化与 near 变化复用当前结果；过滤条件变化后需重新点击读取。"),
      p("没有提供 CRAM、GTF 注释、多样本比较、fragment/UMI 计数、链特异分母或正式 PSI。"))
  )
)

server <- function(input, output, session) {
  current_config <- reactive({
    list(demo = identical(input$source, "demo"),
      bam = if (identical(input$source, "demo")) "" else input$bam,
      chrom = trimws(input$chrom), start1 = input$roi_start, end1 = input$roi_end,
      mapq = input$mapq, baseq = input$baseq, anchor = input$anchor,
      min_intron = input$min_intron, max_intron = input$max_intron,
      strand_mode = input$strand_mode, exclude_duplicates = input$exclude_dup, nh1_only = input$nh1)
  })
  launched_config <- reactiveVal(NULL)
  task <- ExtendedTask$new(function(cfg, source_file) {
    promises::future_promise({
      e <- new.env(parent = globalenv())
      sys.source(source_file, envir = e)
      e$analyze_bam(cfg)
    }, seed = TRUE)
  }) |> bind_task_button("run")
  observeEvent(input$run, {
    cfg <- current_config()
    if (!cfg$demo && (is.null(cfg$bam) || !cfg$bam %in% unname(bams))) {
      showNotification("指定目录没有可用 BAM；请设置 REGTOOLS_BAM_DIR 后重新启动。", type = "error")
      bslib::update_task_button("run", state = "ready", session = session)
      return()
    }
    launched_config(cfg)
    task$invoke(cfg, backend_file)
  })
  result <- reactive(task$result())
  metrics <- reactive({
    r <- result()
    out <- tryCatch(backend$summarize_target(r, input$target_start, input$target_end, input$delta),
                    error = function(e) e)
    validate(need(!inherits(out, "error"), if (inherits(out, "error")) conditionMessage(out) else ""))
    out
  })
  output$status <- renderUI({
    s <- task$status()
    if (s == "initial") return(p("尚未读取。先点击“读取区间”；默认模拟数据的预期结果：4 reads，exact=2，near=1，每100 reads=50。"))
    if (s == "running") return(p("正在计算已提交的参数快照；完成前不展示本次未核验的计数。"))
    if (s == "error") {
      err <- tryCatch(task$result(), error = function(e) conditionMessage(e))
      return(div(class = "alert alert-danger", "计算未完成：", as.character(err)))
    }
    r <- task$result()
    stale <- !identical(current_config(), launched_config())
    tagList(
      div(class = if (stale) "alert alert-warning" else "alert alert-success",
          if (stale) "读取参数已修改：下方仍是上一次已完成结果。点击读取后应用新参数。" else "计算完成，原生 RegTools score 与逐 read 计数核验通过。"),
      if (length(r$warnings)) div(class = "alert alert-warning", paste(r$warnings, collapse = "；")))
  })
  output$metrics <- renderUI({
    m <- metrics()
    num <- function(x, digits = 0L) if (is.na(x)) "NA" else format(round(x, digits), big.mark = ",", trim = TRUE)
    layout_column_wrap(width = 1/3,
      value_box("Junction reads（区间内，并集）", num(m$junction_reads)),
      value_box("Read depth（区间均值）", num(m$mean_read_depth, 3)),
      value_box("Junction exact count", num(m$junction_exact_count)),
      value_box("Junction near count（不含 exact）", num(m$junction_near_count)),
      value_box("Junction per 100 reads", num(m$junction_per_100_reads, 2)))
  })
  output$denominator_note <- renderText({
    m <- metrics()
    paste0("固定分母：", m$denominator_reads, " 条 primary alignment；区间 ",
      m$chrom, ":", m$denominator_start1, "-", m$denominator_end1,
      "。每 100 reads 使用 exact count，不使用 near count 或平均 depth。")
  })
  output$display_control <- renderUI({
    r <- result()
    sliderInput("display", "显示窗口（不改变统计区间）", min = r$config$start1,
                max = r$config$end1, value = c(r$config$start1, r$config$end1), step = 1, sep = "")
  })
  view_range <- reactive({
    r <- result(); x <- input$display
    if (is.null(x) || length(x) != 2L || x[[1L]] < r$config$start1 || x[[2L]] > r$config$end1)
      c(r$config$start1, r$config$end1) else x
  })
  output$depth_plot <- renderPlot({
    r <- result(); w <- view_range()
    d <- r$depth[r$depth$pos1 >= w[[1L]] & r$depth$pos1 <= w[[2L]], , drop = FALSE]
    plot(d$pos1, d$depth, type = "l", xlab = paste0(r$config$chrom, " · 1-based position"),
         ylab = "Read depth", main = "Aligned-base coverage", ylim = c(0, max(1, d$depth)))
    abline(v = c(input$target_start, input$target_end), lty = 2)
  })
  junction_rows <- reactive({
    r <- result(); j <- r$junctions
    j <- j[order(-j$score, j$intron_start1), , drop = FALSE]
    exact <- j$intron_start1 == input$target_start & j$intron_end1 == input$target_end
    near <- abs(j$intron_start1 - input$target_start) <= input$delta &
      abs(j$intron_end1 - input$target_end) <= input$delta & !exact
    j$target_class <- ifelse(exact, "exact", ifelse(near, "near", "other"))
    j
  })
  output$junction_plot <- renderPlot({
    w <- view_range(); j <- junction_rows()
    j <- head(j[j$intron_start1 >= w[[1L]] & j$intron_end1 <= w[[2L]], , drop = FALSE], 30L)
    plot(NA, xlim = w, ylim = c(0, 1.25), xlab = "Intron positions (1-based)", ylab = "",
         yaxt = "n", main = "Native RegTools junction support")
    if (!nrow(j)) { text(mean(w), .5, "No junctions fully inside this display window"); return(invisible(NULL)) }
    for (i in seq_len(nrow(j))) {
      t <- seq(0, 1, length.out = 80)
      height <- .2 + .75 * j$score[[i]] / max(j$score)
      x <- j$intron_start1[[i]] + t * (j$intron_end1[[i]] - j$intron_start1[[i]])
      y <- 4 * height * t * (1 - t)
      lines(x, y, lwd = if (j$target_class[[i]] == "exact") 3 else 1,
            lty = if (j$target_class[[i]] == "near") 2 else 1)
      text(mean(range(x)), height + .035, labels = paste0(j$score[[i]], " ", j$strand[[i]]), cex = .8)
    }
  })
  output$junction_table <- DT::renderDT({
    j <- junction_rows()
    j[, c("chrom", "intron_start1", "intron_end1", "strand", "score", "target_class",
           "intron_start0", "intron_end0", "blockSizes"), drop = FALSE]
  }, rownames = FALSE, selection = "single", options = list(pageLength = 10, scrollX = TRUE))
  observeEvent(input$junction_table_rows_selected, {
    i <- input$junction_table_rows_selected
    j <- junction_rows()
    if (length(i) == 1L && i <= nrow(j)) {
      updateNumericInput(session, "target_start", value = j$intron_start1[[i]])
      updateNumericInput(session, "target_end", value = j$intron_end1[[i]])
    }
  })
  all_evidence <- reactive({
    r <- result(); e <- r$events
    extra <- setdiff(names(r$reads), c("read_id", "chrom", "strand"))
    ans <- cbind(e, r$reads[match(e$read_id, r$reads$read_id), extra, drop = FALSE])
    rownames(ans) <- NULL
    ans
  })
  output$evidence_table <- DT::renderDT({
    e <- all_evidence()
    if (!input$evidence_all) e <- e[
      abs(e$intron_start1 - input$target_start) <= input$delta &
      abs(e$intron_end1 - input$target_end) <= input$delta, , drop = FALSE]
    e[, setdiff(names(e), "key"), drop = FALSE]
  }, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE))
  output$audit <- renderPrint({
    r <- result()
    cat("RegTools score / CIGAR evidence:", if (r$native_audit_passed) "PASS" else "FAIL", "\n")
    cat("候选 alignment:", r$candidate_alignments, "\n")
    cat("仅 reference span 重叠、无 M/= /X 覆盖而排除:", r$span_only_excluded, "\n")
    cat("保留 alignment:", nrow(r$reads), "\n")
    cat("完成时间:", r$completed_at, "\n\n")
    print(metrics())
  })
  output$download_summary <- downloadHandler("junction_metrics.csv", function(file)
    utils::write.csv(metrics(), file, row.names = FALSE, na = "NA", fileEncoding = "UTF-8"))
  output$download_junctions <- downloadHandler("junctions.csv", function(file)
    utils::write.csv(junction_rows(), file, row.names = FALSE, fileEncoding = "UTF-8"))
  output$download_depth <- downloadHandler("depth.tsv", function(file)
    utils::write.table(result()$depth, file, sep = "\t", quote = FALSE, row.names = FALSE))
  output$download_evidence <- downloadHandler("junction_read_evidence.csv", function(file)
    utils::write.csv(all_evidence(), file, row.names = FALSE, na = "NA", fileEncoding = "UTF-8"))
  output$download_metadata <- downloadHandler("parameters_and_audit.json", function(file) {
    r <- result()
    x <- list(config = r$config, target_metrics = metrics(), versions = r$versions,
      source_metadata = r$source_metadata, completed_at = r$completed_at,
      native_audit_passed = r$native_audit_passed, warnings = r$warnings, commands = r$log,
      definitions = list(count_unit = "primary alignment records, paired ends separate",
        denominator = "retained reads with M/= /X blocks overlapping fixed analysis ROI",
        near = "both boundaries within delta; exclude exact events; union read IDs",
        per100 = "100 * exact / denominator; NA when denominator is zero",
        depth = "samtools depth; zero-inclusive mean; BQ affects depth only; pair overlaps counted twice"))
    jsonlite::write_json(x, file, auto_unbox = TRUE, pretty = TRUE, na = "null")
  })
}
shinyApp(ui, server)
