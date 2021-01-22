# Run from the main folder as  ./Scripts/download_ncbi.sh $ACC $col 
# or ./Scripts/download_ncbi.sh $ACC      
# ($ACC is the accession number for example CP001312.1)
# ($col is the folder where you want the data to be saved, default:QUERY)

ACC=$1
col="QUERY"

if [ -n "$2" ]; then
col=$2
fi

mkdir ${col}/${ACC}
wget -O ${col}/${ACC}/${ACC}.fna "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=fasta&id=${ACC}"
wget -O ${col}/${ACC}/${ACC}.gff "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=${ACC}"

# If you don't have wget you can use curl
# curl -o ${col}/${ACC}/${ACC}.fna "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=fasta&id=${ACC}"
# curl -o ${col}/${ACC}/${ACC}.gff "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=${ACC}"
