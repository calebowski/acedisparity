args <- commandArgs(trailingOnly = TRUE)
replicate_id <- as.numeric(args[1])
tree_size <- args[2]
job_id <- as.numeric(args[3])


base_path <- paste0("/mnt/parscratch/users/bip24cns/acedisparity/randomExtinction/continuous/", tree_size, "/")
# job_id <- 8469104
write.path <- function(subfolder, filename) {
  return(paste0(base_path, subfolder, "/", job_id, "_", sprintf(filename, replicate_id)))
}


ord_sample <- list(bm = readRDS(write.path("ord/temp", paste0("bm_ord_sample_%03d.rds"))),
                   bm_t = readRDS(write.path("ord/temp", paste0("bm_t_ord_sample_%03d.rds"))),
                   ou_w = readRDS(write.path("ord/temp", paste0("ou_w_ord_sample_%03d.rds"))),
                   ou_st = readRDS(write.path("ord/temp", paste0("ou_st_ord_sample_%03d.rds"))),
                   ou_sh = readRDS(write.path("ord/temp", paste0("ou_sh_ord_sample_%03d.rds")))
)


ord_point <- list(bm = readRDS(write.path("ord/temp", paste0("bm_ord_point_%03d.rds"))),
                   bm_t = readRDS(write.path("ord/temp", paste0("bm_t_ord_point_%03d.rds"))),
                   ou_w = readRDS(write.path("ord/temp", paste0("ou_w_ord_point_%03d.rds"))),
                   ou_st = readRDS(write.path("ord/temp", paste0("ou_st_ord_point_%03d.rds"))),
                   ou_sh = readRDS(write.path("ord/temp", paste0("ou_sh_ord_point_%03d.rds")))
)

ord_no_ace <- list(bm = readRDS(write.path("ord/temp", paste0("bm_ord_no_ace_%03d.rds"))),
                   bm_t = readRDS(write.path("ord/temp", paste0("bm_t_ord_no_ace_%03d.rds"))),
                   ou_w = readRDS(write.path("ord/temp", paste0("ou_w_ord_no_ace_%03d.rds"))),
                   ou_st = readRDS(write.path("ord/temp", paste0("ou_st_ord_no_ace_%03d.rds"))),
                   ou_sh = readRDS(write.path("ord/temp", paste0("ou_sh_ord_no_ace_%03d.rds")))
)


ord_true <- list(bm = readRDS(write.path("ord/temp", paste0("bm_ord_true_%03d.rds"))),
                   bm_t = readRDS(write.path("ord/temp", paste0("bm_t_ord_true_%03d.rds"))),
                   ou_w = readRDS(write.path("ord/temp", paste0("ou_w_ord_true_%03d.rds"))),
                   ou_st = readRDS(write.path("ord/temp", paste0("ou_st_ord_true_%03d.rds"))),
                   ou_sh = readRDS(write.path("ord/temp", paste0("ou_sh_ord_true_%03d.rds")))
)

post_ord_point <- list(bm = readRDS(write.path("ord/temp", paste0("bm_post_ord_point_%03d.rds"))),
                   bm_t = readRDS(write.path("ord/temp", paste0("bm_t_post_ord_point_%03d.rds"))),
                   ou_w = readRDS(write.path("ord/temp", paste0("ou_w_post_ord_point_%03d.rds"))),
                   ou_st = readRDS(write.path("ord/temp", paste0("ou_st_post_ord_point_%03d.rds"))),
                   ou_sh = readRDS(write.path("ord/temp", paste0("ou_sh_post_ord_point_%03d.rds")))
)

post_ord_sample <- list(bm = readRDS(write.path("ord/temp", paste0("bm_post_ord_sample_%03d.rds"))),
                   bm_t = readRDS(write.path("ord/temp", paste0("bm_t_post_ord_sample_%03d.rds"))),
                   ou_w = readRDS(write.path("ord/temp", paste0("ou_w_post_ord_sample_%03d.rds"))),
                   ou_st = readRDS(write.path("ord/temp", paste0("ou_st_post_ord_sample_%03d.rds"))),
                   ou_sh = readRDS(write.path("ord/temp", paste0("ou_sh_post_ord_sample_%03d.rds")))
)

saveRDS(ord_sample, write.path("ord", "ord_sample_%03d.rds"))
saveRDS(ord_point, write.path("ord", "ord_point_%03d.rds"))
saveRDS(ord_no_ace, write.path("ord", "ord_no_ace_%03d.rds"))
saveRDS(ord_true, write.path("ord", "ord_true_%03d.rds"))
saveRDS(post_ord_point, write.path("ord", "post_ord_point_%03d.rds"))
saveRDS(post_ord_sample,  write.path("ord", "post_ord_sample_%03d.rds"))