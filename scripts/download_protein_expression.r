##download_protein_expression.r
##2015-01-13 dmontaner@cipf.es
##2015-04-23 julenmendieta92@gmail.com
##Collecting data from TCGA

## The scripts uses TCGA DCC Web Services to find out all the PROTEIN EXPRESSION data.

## Batch effect is very strong.
## Do just comparisons whiting each disease.

## Find the platform here:
## http://tcga-data.nci.nih.gov/tcgadccws/GetHTML?query=Platform

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

#setwd (file.path (.job$dir$raw, "protein_exp"))
setwd (file.path ("/home/jmendieta/Documents/tcga_download/Download/expre"))

################################################################################

## DOWNLOAD DATA

base   <- "https://tcga-data.nci.nih.gov"
my.url <- "http://tcga-data.nci.nih.gov/tcgadccws/GetXML?query=Archive[Platform[@name=MDA_RPPA_Core]][ArchiveType[@type=Level_3]][@deployStatus=Available][@isLatest=1]"

char <- getURL (my.url)
char

mix <- xmlInternalTreeParse (char)
class (mix)
mix <- xmlChildren (xmlChildren (mix)[["httpQuery"]])[["queryResponse"]]
class (mix)
mix

ns <- getNodeSet (mix, path = "class")
class (ns)
ns

datos <- xmlToDataFrame (ns, stringsAsFactors = FALSE)
dim (datos)
colnames (datos) <- paste0 ("V", 1:ncol (datos))
colnames (datos)[c (1, 3, 7)] <- c("date", "location", "name")
datos[1:3, 1:8]

datos[,"fecha"] <- as.Date (datos[,"date"], "%m-%d-%Y")
datos[,c ("date", "fecha")]

datos[,"url"] <- paste0 (base, datos[,3])
##datos[,"gzfile"] <- basename (datos[,3])
##datos[,"gzfolder"]  <- sub (".tar.gz", "", datos[,"gzfile"])

rownames (datos) <- datos[,"name"]

datos[1:3,]

system.time ({
    for (url in datos[,"url"]) {
        filename <- basename (url)
        download.file (url = url, destfile = filename, method = "curl")
        untar (tarfile = filename, compressed = TRUE)    
    }
})

################################################################################

## SOME CHECKS
length (dir (recursive= TRUE))

ficheros <- dir (pattern = "protein_expression", recursive= TRUE)
length (ficheros)
ficheros[1:3]

li <- strsplit (ficheros, ".protein_expression.Level_3.")
table (sapply (li, length))  ##OK all 2
uuid <- sapply (li, function (x) x[2])
uuid <- sub (".txt", "", uuid)
length (uuid)

folder <- dirname (ficheros)

fi <- sub ("mdanderson.org_", "", ficheros)
li <- strsplit (fi, "\\.")
table (sapply (li, length))  ##OK all the same
tag <- sapply (li, function (x) x[1])
length (tag)

finfo <- as.data.frame (list (fichero = ficheros, uuid = uuid, folder = folder, tag = tag), stringsAsFactors= FALSE)
class (finfo)
dim (finfo)
sapply (finfo, class)
table (finfo[,"tag"])

finfo[,"uuid.infile"] <- NA
system.time ({
    for (i in 1:nrow (finfo)) {
        finfo[i, "uuid.infile"] <- readLines (con = finfo[i, "fichero"], n = 1)
    }
})
finfo[,"uuid.infile"] <- sub ("Sample REF\t", "", finfo[,"uuid.infile"])
table (finfo[,"uuid"] == finfo[,"uuid.infile"], exclude = NULL)  ## OK all the same

if (any (finfo[,"uuid"] != finfo[,"uuid.infile"])) stop ("UUIDs do not match")

summary (finfo)

########################################

## Some duplicated UUIDS ????
dups <- duplicated (finfo[,"uuid"])
table (dups)

duplicados <- finfo[dups, "uuid"]
finfo[,"is.dup"] <- finfo[,"uuid"] %in% duplicados
table (finfo[,"is.dup"], exclude = NULL)

table (table (finfo[,"uuid"], exclude = NULL), exclude = NULL)
table (finfo[,"tag"], finfo[,"is.dup"])

##sort by date and SERIAL INDEX
table (finfo$folder %in% rownames (datos))  ## OK
finfo[,"fecha"] <- datos[finfo$folder, "fecha"]
orden <- order (finfo[,"fecha"], finfo[,"is.dup"], finfo[,"uuid"], finfo[,"folder"], decreasing = TRUE)
finfo <- finfo[orden,]
finfo[finfo$is.dup,][1:10,]

malos <- finfo[finfo$is.dup,]
orden <- order (malos[,"uuid"])
malos <- malos[orden,]
table (malos[,"tag"])
malos[1:10,]
tail (malos)

dup <- duplicated (finfo[,"uuid"])
table (dup)
finfo <- finfo[!dup,]
table (duplicated (finfo[,"uuid"]))
     
rownames (finfo) <- finfo[,"uuid"]

################################################################################

## READ DATA
datos.li <- list ()
system.time ({
    for (i in 1:nrow (finfo)) {
        datos.li[[finfo[i, "uuid"]]] <-
            read.table (file = finfo[i, "fichero"], header = TRUE, sep = "\t", quote = "", skip = 1, colClasses = c ("character", "numeric"), row.names = 1)
    }
})
length (datos.li)
table (names (datos.li) == rownames (finfo))

unique (sapply (datos.li, colnames)) ## OK all columns the same


## MAKE MATRIX
id.list <- lapply (datos.li, rownames)
summary (sapply (id.list, length))
summary (sapply (datos.li, nrow))

ids <- sort (unique (unlist (id.list)))
length (ids)
ids

prot.exp <- matrix (NA, nrow = length (ids), ncol = length (datos.li))
dim (prot.exp)
rownames (prot.exp) <- ids
colnames (prot.exp) <- names (datos.li)
prot.exp[1:3, 1:4]
system.time ({
    for (uid in names (datos.li)) {
        mat <- datos.li[[uid]]
        prot.exp[,uid] <- mat[ids,]
    }
})

prot.exp[1:3, 1:4]

table (colnames (prot.exp) == rownames (finfo))

################################################################################

## ## Some Plots
## table (finfo$tag)

## boxplot (prot.exp[,c(which (finfo$tag == "ACC"), which (finfo$tag == "BLCA"))])
## abline (h = 0, col = "blue")
## boxplot (prot.exp[,c(which (finfo$tag == "KIRP"), which (finfo$tag == "KIRC"))])
## abline (h = 0, col = "blue")

## orden <- order (finfo$tag)
## boxplot (prot.exp[,orden])
## abline (h = 0, col = "blue")

## ran <- sample (colnames (prot.exp), size = 100)
## boxplot (prot.exp[,ran])
## abline (h = 0, col = "blue")

################################################################################

### SAVE
prot.exp.info <- finfo
table (colnames (prot.exp) == rownames (prot.exp.info))
#save (list = c("prot.exp", "prot.exp.info"), file = file.path (.job$dir$proces, "prot_exp.RData"))
save (list = c("prot.exp", "prot.exp.info"), file = file.path ("/home/jmendieta/Documents/tcga_download/Download/prot_exp.RData"))


###EXIT
warnings ()
sessionInfo ()
q ("no")
