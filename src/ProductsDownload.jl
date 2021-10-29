module ProductsDownload
"""
Module for the download of satellitar data of the Copernicus project
"""



using CSV
using Dates
using DataFrames
using Downloads
using HTTP
using NCDatasets
using ZipFile


include("ProductsData.jl")


test = [ "C:\\Users\\DAVIDE-FAVARO\\Desktop\\XML\\1.xml",
         "C:\\Users\\DAVIDE-FAVARO\\Desktop\\XML\\2.xml",
         "C:\\Users\\Lenovo\\Desktop\\XML\\1.xml",
         "C:\\Users\\Lenovo\\Desktop\\XML\\2.xml",
         "C:\\Users\\Lenovo\\Desktop\\XML\\Prod_Test.xml" ]

out = [ "D:\\Vario\\Stage",
        "C:\\Users\\Lenovo\\Desktop\\XML",
        "C:\\Users\\Lenovo\\Desktop\\XML\\Test",
        "C:\\Users\\DAVIDE-FAVARO\\Desktop",
        "C:\\Users\\DAVIDE-FAVARO\\Desktop\\XML" ]




# https://scihub.copernicus.eu/dhus/odata/v1
# https://scihub.copernicus.eu/dhus/search?q=*&rows=25
# https://scihub.copernicus.eu/dhus/odata/v1/Products('2b17b57d-fff4-4645-b539-91f305c27c69')/$value

# S3B_SL_2_LST____20210915T101845_20210915T102145_20210915T124701_0179_057_065_2160_LN2_O_NR_004  Nome di un file, tentare di usarlo in Product(...) ritorna not found
# 1def7d25-ecd6-49ef-8e3a-243f35857951 uuid file 240 Mb
# b0526ddf-c3f3-4065-8dc5-1891ffc8e326 uuid file 165 Mb

# https://scihub.copernicus.eu/dhus/odata/v1/Products('2b17b57d-fff4-4645-b539-91f305c27c69')/$value
# https://scihub.copernicus.eu/dhus/odata/v1/Products?filter=Name eq 'S3B_SL_2_LST____20210915T101845_20210915T102145_20210915T124701_0179_057_065_2160_LN2_O_NR_004'



"""
    getData( fileId::AbstractString, targetDirectory::AbstractString, authToken::AbstractString )

Download the archive containing the data identified by "fileId", into the chosen directory using the user authentication token
"""
function getData( fileId::AbstractString, targetDirectory::AbstractString, authToken::AbstractString )                 

    Downloads.download( "https://scihub.copernicus.eu/dhus/odata/v1/Products('$fileId')/\$value",
                        targetDirectory*"\\$(fileId[1:10]).zip", #Destinazione
                        headers = [ "Authorization" => "Basic $authToken" ] #=,
                        progress = ( total, now ) -> println("$(now/total*100)% ( $now / $total )"), #Funzione che permette di controllare lo stato del download
                        verbose = true =# ) #Più info) #info di autorizzazioone

    
end

#   getData( "b57f225e-d288-4e4a-bb35-2a7eb75d60e4", out[3], authenticate("davidefavaro", "Tirocinio")  )



"""
    unzipIn( dir::AbstractString )

Unzip a ".zip" file stored in "dir", creating a new directory within "dir" containing the file contents 
"""
function unzip( dir::AbstractString )
    zipPath = "$dir\\$(readdir(dir)[1])"
    println(zipPath)
    zip = ZipFile.Reader(zipPath)
    new = mkdir(zipPath[1:end-4])
    newdir = zipPath[1:end-4]*"\\"
    start = length(zip.files[1].name) + 1
    for i in 2:length(zip.files)
        write( newdir*zip.files[i].name[start:end], read(zip.files[i]) )
    end
    return new
end

#   dir = unzip( out[3] )