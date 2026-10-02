source("backend.R"); source("variant_backend.R"); source("variant_jobs.R")
expect_error <- function(expr, pattern=NULL) {
  x<-tryCatch({force(expr);NULL},error=function(e)conditionMessage(e));stopifnot(!is.null(x))
  if(!is.null(pattern))stopifnot(grepl(pattern,x,fixed=TRUE))
}
local({
  root<-tempfile("regshiny_jobs_");dir.create(root,mode="0700");on.exit(unlink(root,recursive=TRUE),add=TRUE)
  bin<-file.path(root,"bin");dir.create(bin);mode<-file.path(root,"scheduler_mode");writeLines("live",mode)
  script<-function(name,body){p<-file.path(bin,name);writeLines(c("#!/bin/sh",body),p);Sys.chmod(p,"0700")}
  script("qsub","printf '12345\\n'")
  script("qstat",paste0("if [ \"$(cat ",shQuote(mode),")\" = live ]; then exit 0; else exit 1; fi"))
  script("qacct",c(paste0("case \"$(cat ",shQuote(mode),")\" in"),
    "success) printf 'jobnumber 12345\\nfailed 0\\nexit_status 0\\n';;",
    "failure) printf 'jobnumber 12345\\nfailed 0\\nexit_status 1\\n';;","*) exit 1;; esac"))
  oldpath<-Sys.getenv("PATH");Sys.setenv(PATH=paste(bin,oldpath,sep=.Platform$path.sep));on.exit(Sys.setenv(PATH=oldpath),add=TRUE)
  oldroot<-Sys.getenv("REGSHINY_JOB_ROOT");Sys.setenv(REGSHINY_JOB_ROOT=root);on.exit(Sys.setenv(REGSHINY_JOB_ROOT=oldroot),add=TRUE)
  uidroot<-variant_user_job_root(root);stopifnot(as.integer(file.info(uidroot)$mode)%%512L==448L)
  outside<-file.path(root,"outside");dir.create(outside,mode="0700");link<-file.path(uidroot,"run-link");file.symlink(outside,link)
  expect_error(validate_variant_job_dir(link,root),"symbolic");expect_error(validate_variant_job_dir(outside,root),"own Reg_Shiny")
  lock<-job_acquire_submit_lock(root);expect_error(job_acquire_submit_lock(root),"already in progress");unlink(lock,recursive=TRUE)
  vcf<-file.path(root,"synthetic.vcf.gz");writeLines("synthetic",vcf);writeLines("index",paste0(vcf,".tbi"))
  bams<-file.path(root,paste0("sample",seq_len(36),".bam"));for(p in bams){writeLines("synthetic",p);writeLines("index",paste0(p,".bai"))}
  registry<-file.path(root,"vcf.tsv");manifest<-file.path(root,"samples.tsv")
  write.table(data.frame(chrom="chrDemo",vcf=vcf,build="SYNTHETIC"),registry,sep="\t",row.names=FALSE,quote=FALSE)
  write.table(data.frame(vcf_sample=paste0("S",seq_len(36)),bam=bams,build="SYNTHETIC"),manifest,sep="\t",row.names=FALSE,quote=FALSE)
  resources<-load_variant_resources(registry,manifest)
  v<-data.frame(record_id="SYNTHETIC:chrDemo:190:A:G",chrom="chrDemo",pos1=190L,ref="A",alt="G",build="SYNTHETIC",vcf=vcf)
  calls<-data.frame(vcf_sample=paste0("S",seq_len(36)),raw_gt=rep(c("0/0","0/1","1/1"),each=12),genotype=rep(c("0/0","0/1","1/1"),each=12),bam=bams,linked=TRUE,call_status="CALLED")
  originals<-list(lookup_variants=lookup_variants,variant_vcf_info=variant_vcf_info,variant_get_genotypes=variant_get_genotypes)
  assign("lookup_variants",function(resources,query)v,.GlobalEnv)
  assign("variant_vcf_info",function(path,chrom)list(state=variant_file_state(c(path,paste0(path,".tbi")))),.GlobalEnv)
  assign("variant_get_genotypes",function(resources,variant)calls,.GlobalEnv)
  on.exit(for(n in names(originals))assign(n,originals[[n]],.GlobalEnv),add=TRUE)
  cfg<-default_demo_config();cfg$demo<-FALSE
  receipt<-submit_variant_job(resources,v,cfg,normalizePath("backend.R"),normalizePath("variant_backend.R"),normalizePath("variant_jobs.R"),root,ui_controls=list(sampling_mode="all"))
  request<-job_request(receipt$job_dir);stopifnot(receipt$job_id=="12345",nrow(request$selected_plan)==36L,length(request$bam_states)==36L,length(request$source_hashes)==3L)
  expect_error(submit_variant_job(resources,v,cfg,"backend.R","variant_backend.R","variant_jobs.R",root),"already active")
  expect_error(load_variant_job_result(receipt$job_dir,root),"not yet verified")
  saveRDS(list(synthetic=TRUE),file.path(receipt$job_dir,"result.rds"))
  job_atomic_json(list(status="APPLICATION_COMPLETE",job_id="12345",failed_samples=0,result_md5=unname(tools::md5sum(file.path(receipt$job_dir,"result.rds")))),file.path(receipt$job_dir,"application_receipt.json"))
  stopifnot(!read_variant_job(receipt$job_dir,root)$result_ready)
  drive<-function(expected){for(i in seq_len(100L)){
    for(k in ls(.variant_job_probes)){p<-get(k,.variant_job_probes);p$last_checked<-0;assign(k,p,.variant_job_probes)}
    z<-read_variant_job(receipt$job_dir,root);if(z$status==expected)return(z);Sys.sleep(.025)
  };stop("Mock state not reached: ",expected)}
  writeLines("transient",mode);drive("ACCOUNTING_PENDING");writeLines("live",mode);drive("QUEUED")
  stopifnot(!read_variant_job(receipt$job_dir,root)$result_ready)
  writeLines("success",mode);done<-drive("COMPLETE");stopifnot(done$result_ready,done$scheduler_terminal,load_variant_job_result(receipt$job_dir,root)$synthetic)
  saveRDS(list(synthetic=FALSE),file.path(receipt$job_dir,"result.rds"));expect_error(load_variant_job_result(receipt$job_dir,root),"not yet verified")
  stopifnot(is.null(job_parse_accounting(c("jobnumber 999","failed 0","exit_status 0"),"12345")))
  s<-job_read_json(file.path(receipt$job_dir,"state.json"));s$status<-"FAILED";s$result_ready<-FALSE;s$accounting<-NULL;job_atomic_json(s,file.path(receipt$job_dir,"state.json"))
  writeLines("live",mode);expect_error(resume_variant_job(receipt$job_dir,root),"terminal scheduler")
  writeLines("failure",mode);failed<-drive("INTERRUPTED");stopifnot(failed$scheduler_terminal,!failed$result_ready)
  resumed<-resume_variant_job(receipt$job_dir,root);stopifnot(resumed$job_id=="12345",file.exists(file.path(receipt$job_dir,"submission-2.json")))
  script("qsub","exit 1");rejected<-file.path(uidroot,"run-rejected");dir.create(rejected,mode="0700");saveRDS(request,file.path(rejected,"request.rds"))
  expect_error(job_submit_saved(rejected,root));stopifnot(job_read_json(file.path(rejected,"state.json"))$status=="SUBMISSION_REJECTED")
  for(k in ls(.variant_job_probes)){p<-get(k,.variant_job_probes);if(!is.null(p$process)&&p$process$is_alive())p$process$kill_tree()}
})
cat("Persistent job isolation, source freeze, accounting gate and resume tests: PASS\n")
