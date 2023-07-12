
library(bnlearn)

args <- commandArgs(trailingOnly = TRUE)

input_filename = args[1]
algo = args[2]

output_filename = paste(algo, 'boot',input_filename, sep = '_')
output_RData = paste(output_filename, '_boot.RData', sep = '')

mat_old <- read.table(input_filename, header = TRUE)

boot <- boot.strength(mat_old, R = 500, algorithm = algo)

boot[(boot$strength > 0.1) & (boot$direction >= 0.5), ]

boot_tmp <- averaged.network(boot, threshold = 0.1)
save.image(output_RData)
write.dot(boot_tmp, file=output_filename)


