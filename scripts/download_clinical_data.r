##download_clinical_data.r
##2015-01-13 dmontaner@cipf.es
##2015-04-23 julenmendieta92@gmail.com
##Collecting data from TCGA

## The scripts uses TCGA DCC Web Services to find out all CLINICAL data.

## Find the platform here:
## http://tcga-data.nci.nih.gov/tcgadccws/GetHTML?query=Platform
## "bio" seems to be the appropriated one

## grep -B 1 -A 2 E6B279BC-97A1-4B9A-9632-06326C418E1C nationwidechildrens.org_biospecimen_sample_coad.txt 

date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ##"R version 3.2.0 (2015-04-16)"
library (RCurl); packageDescription ("RCurl", fields = "Version") #"1.95-4.5"
library (XML); packageDescription ("XML", fields = "Version") #"3.98-1.1"
#help (package = XML)
#help (package = RCurl)

try (source (".job.r")); try (.job)

options (width = 170)
#options (width = 1000)

#setwd (file.path (.job$dir$raw, "clinical"))
setwd (file.path ("/home/jmendieta/Documents/tcga_download/Download/Clinicos"))


################################################################################

## DOWNLOAD DATA
base   <- "https://tcga-data.nci.nih.gov"
my.url <- "http://tcga-data.nci.nih.gov/tcgadccws/GetXML?query=Archive[Platform[@name=bio]][ArchiveType[@type=Level_2]][@deployStatus=Available][@isLatest=1]"

char <- getURL (my.url)
char

mix <- xmlInternalTreeParse (char)
class (mix)
mix <- xmlChildren (xmlChildren (mix)[["httpQuery"]])[["queryResponse"]]
class (mix)
mix

ns <- getNodeSet (mix, path = "class")
class (ns)
length (ns)
ns

datos <- xmlToDataFrame (ns, stringsAsFactors = FALSE)
dim (datos)
datos[1:3,]
colnames (datos) <- paste0 ("V", 1:ncol (datos))
colnames (datos)[c (1, 2, 3, 7)] <- c("date", "baseName", "location", "name")
datos[1:3,]

datos[,"fecha"] <- as.Date (datos[,"date"], "%m-%d-%Y")
datos[,c ("date", "fecha")]

datos[,"url"] <- paste0 (base, datos[,3])

datos[,"tag"] <- sapply (strsplit (datos[,"baseName"], split = "_"), function (x) x[2])
revdups <- table (duplicated (datos[,"tag"]))

##Si hay duplicados los eliminamos para que el siguiente paso no de error
if (revdups[2] >= 1) {dup <- duplicated (datos[,"tag"]) ; datos <- datos[!dup,]} 

rownames (datos) <- datos[,"tag"]
datos[1:3,]

system.time ({
    for (url in datos[,"url"]) {
        filename <- basename (url)
        download.file (url = url, destfile = filename, method = "curl")
        untar (tarfile = filename, compressed = TRUE)    
    }
})

################################################################################

## FIND UUID
ficheros <- dir (pattern = "biospecimen_shipment_portion", recursive= TRUE)
length (ficheros)

datos <- NULL
for (fi in ficheros) {
    print (fi)
    nombres <- unlist (strsplit (readLines (fi)[1], split = "\t"))
    #print (nombres)
    datos0 <- read.table (fi, header = FALSE, sep = "\t", quote = "", col.names = nombres, as.is = TRUE, na.strings = "[Not Available]", skip = 2)
    datos0[,"file"] <- fi
    print (dim (datos0))
    try (datos <- rbind (datos, datos0))
}
dim (datos)

datos[,"tag"] <- sapply (strsplit (sapply (strsplit (datos[,"file"], split = "_"), function (x) x[2]), split = "\\."), function (x) x[1])
table (datos[,"tag"], exclude = NULL)

if (any (duplicated (datos[,"bcr_shipment_portion_uuid"]))) stop ("DUPLICATED") ## OK NO DUPS

##protein data
#load (file.path (.job$dir$proces, "prot_exp.RData"))
load (file.path ("/home/jmendieta/Documents/tcga_download/Download/prot_exp.RData"))

dim (prot.exp)

prot.ids <- colnames (prot.exp)
prot.ids[1:3]

table (datos[,"bcr_shipment_portion_uuid"] %in% prot.ids)
table (prot.ids %in% datos[,"bcr_shipment_portion_uuid"]) #deberian esta todos ???
table (tolower (prot.ids) %in% tolower (datos[,"bcr_shipment_portion_uuid"])) #deberian esta todos ???
setdiff (tolower (prot.ids), tolower (datos[,"bcr_shipment_portion_uuid"])) #deberian esta todos ???

malos <- !tolower (prot.ids) %in% tolower (datos[,"bcr_shipment_portion_uuid"])
table (malos)
table (prot.exp.info[malos, "tag"])

datos[,"has.prot"] <- tolower (datos[,"bcr_shipment_portion_uuid"]) %in% tolower (prot.ids)
table (datos[,"has.prot"])
table (datos[,"tag"], datos[,"has.prot"])

################################################################################

## FIND CLASSES
ficheros <- dir (pattern = "biospecimen_sample", recursive = TRUE)
length (ficheros)
ficheros

##apagno un fichero que esta mal
#No se si aqui estara mal tambien, el level a pasado de 2.0.17 a 2.0.25
li <- readLines ("nationwidechildrens.org_COAD.bio.Level_2.0.25.0/nationwidechildrens.org_biospecimen_sample_coad.txt")
li[1:10]
lisp <- strsplit (li, split = "\t")
lon <- sapply (lisp, length)
table (lon)
which (lon == 20)

#Ya no hay este error
#li[376] <- sub ("E6B279BC-97A1-4B9A-9632-06326C418E1C", "[Not Available]\tE6B279BC-97A1-4B9A-9632-06326C418E1C", li[376])
#writeLines (li, con = "nationwidechildrens.org_COAD.bio.Level_2.0.25.0/nationwidechildrens.org_biospecimen_sample_coad.txt")


datos.li <- list ()
for (fi in ficheros) {
    cat ("\n======================== ", fi, " ==============================\n")
    tag <- sub (".txt", "", rev (unlist (strsplit (fi, split = "_")))[1])
    ##
    nombres <- unlist (strsplit (readLines (fi)[1], split = "\t"))
    datos0 <- read.table (fi, header = FALSE, sep = "\t", quote = "", col.names = nombres, as.is = TRUE, na.strings = "[Not Available]", skip = 2)
    datos0[,"tag"] <- tag
    datos.li[[fi]] <- datos0
}


t (sapply (datos.li, dim))
colSums (t (sapply (datos.li, dim)))

unicol <- unique (unlist (sapply (datos.li, colnames)))
datos1 <- NULL
for (fi in ficheros) {
    cat ("\n======================== ", fi, " ==============================\n")
    datos0 <- datos.li[[fi]]
    faltan <- setdiff (unicol, colnames (datos0))
    datos0[,faltan] <- NA
    datos0 <- datos0[,unicol]
    datos1 <- rbind (datos1, datos0)
}
dim (datos1)

sapply (datos1, class)

datos1[1:3,]

table (datos1[,"tag"], exclude = NULL)
table (datos1[,"sample_type"], exclude = NULL)

################################################################################

## COMBINE
dim (datos)
dim (datos1)

datos1[1:3,]
datos[1:3,]

datos1[1:3, c ("bcr_sample_barcode", "bcr_sample_uuid", "sample_type")]
datos [1:3, c ("bcr_sample_barcode", "shipment_portion_bcr_aliquot_barcode", "bcr_shipment_portion_uuid")]

intersect (colnames (datos), colnames (datos1))

table (datos[,"bcr_sample_barcode"] %in% datos1[,"bcr_sample_barcode"])
setdiff (datos[,"bcr_sample_barcode"], datos1[,"bcr_sample_barcode"])
table (datos1[,"bcr_sample_barcode"] %in% datos[,"bcr_sample_barcode"])

combi <- merge (datos, datos1, by = "bcr_sample_barcode", all = TRUE)
dim (combi)
combi[1:3,]

table (combi$has.prot, exclude = NULL)
combi[is.na (combi$has.prot), "has.prot"] <- FALSE

table (tolower (combi[,"tag.x"]) == tolower (combi[,"tag.y"]), exclude = NULL)
igual <- tolower (combi[,"tag.x"]) == tolower (combi[,"tag.y"])
combi[which (!igual), c("tag.x", "tag.y")]


table (combi[combi$has.prot, c("sample_type", "tag.x")])
table (combi[combi$has.prot, c("sample_type", "tag.y")])

###EXIT
warnings ()
sessionInfo ()
q ("no")
