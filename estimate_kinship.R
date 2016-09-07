library(GENESIS)
library(SNPRelate)

# Create GDS file (only do once)
snpgdsBED2GDS(bed.fn='/home/unix/armartin/finrisk/data/popgen_genotypes/imputed_data/Engagex_Egfext_Migraine_FinnishCardio_FTC_geno05_maf05_ld50.bed',
              bim.fn='/home/unix/armartin/finrisk/data/popgen_genotypes/imputed_data/Engagex_Egfext_Migraine_FinnishCardio_FTC_geno05_maf05_ld50.bim',
              fam.fn='/home/unix/armartin/finrisk/data/popgen_genotypes/imputed_data/Engagex_Egfext_Migraine_FinnishCardio_FTC_geno05_maf05_ld50.fam',
              out.gdsfn='/home/unix/armartin/finrisk/data/popgen_genotypes/imputed_data/Engagex_Egfext_Migraine_FinnishCardio_FTC_geno05_maf05_ld50.gds')

geno <- snpgdsOpen('/home/unix/armartin/finrisk/data/popgen_genotypes/imputed_data/Engagex_Egfext_Migraine_FinnishCardio_FTC_geno05_maf05_ld50.gds')

geno <- GdsGenotypeReader(filename = "/home/unix/armartin/finrisk/data/popgen_genotypes/imputed_data/Engagex_Egfext_Migraine_FinnishCardio_FTC_geno05_maf05_ld50.gds")
genoData <- GenotypeData(geno)
mypcair <- pcair(geno)
KINGmat <- snpgdsIBDKING(geno)

iids <- getScanID(geno)

# create matrix of KING estimates
KINGmat <- king2mat(file.kin0 = "/home/unix/armartin/finrisk/data/popgen_genotypes/imputed_data/Engagex_Egfext_Migraine_FinnishCardio_FTC_geno05_maf05_ld50.kin0", 
                    file.kin = "/home/unix/armartin/finrisk/data/popgen_genotypes/imputed_data/Engagex_Egfext_Migraine_FinnishCardio_FTC_geno05_maf05_ld50.kin", 
                    iids = iids)
mypcair <- pcair(genoData = genoData, kinMat = KINGmat, divMat = KINGmat)