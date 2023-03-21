args=commandArgs(TRUE)
ref_seg_BAF = as.numeric(args[1])
ref_seg_LogR = as.numeric(args[2])
ref_seg_nMajor = as.integer(args[3])
ref_seg_nMinor = as.integer(args[4])


gamma_param = 1 # 0.55 for array data

rho = (2*ref_seg_BAF - 1) / (2*ref_seg_BAF - ref_seg_BAF*(ref_seg_nMajor + ref_seg_nMinor) - 1 + ref_seg_nMajor)
psi = (rho*(ref_seg_nMajor + ref_seg_nMinor) + 2 - 2*rho)/(2^(ref_seg_LogR/gamma_param))
psi = (psi-2*(1-rho))/rho

print(paste("rho from reference segment:", rho))
print(paste("psi from reference segment:", psi))

