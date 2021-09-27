module App





using Genie
using SQLite
using CSV
using DataFrames

Genie.newapp_mvc("CopernicusData")


include(joinpath("config", "initializers", "searchlight.jl"))
SQLite.DB("db/data.sqlite")



Genie.newresource("product")





df = CSV.read( "D:\\Documents and Settings\\DAVIDE-FAVARO\\My Documents\\GitHub\\Tirocinio\\Xml di prova\\data.csv", DataFrame )
cols = names(df)
types = eltype.( eachcol(df) )
# SQLite.createtable!( SQLite.DB("db/data_db.sqlite"), "products", Tables.Schema(cols, types) )




io = open( "D:\\Documents and Settings\\DAVIDE-FAVARO\\My Documents\\GitHub\\Tirocinio\\App\\colonne.txt", "w" )
for i in 1:82
    println( io, "column( :$(cols[i]), $(types[i]) )" )
end
close( io )






end # module

#=
column( :uuid, String )
      column( :relpassnumber, Int64 ) #Nullable
      column( :relpassdirection, String ) #Nullable
      column( :filename, String )
      column( :instrumentname, String )
      column( :timeliness, String ) #Nullable
      column( :footprint, String )
      column( :size, String )
      column( :productlevel, String ) #Nullable
      column( :relorbitdir, String ) #Nullable
      column( :ingestiondate, String )
      column( :ecmwf, String ) #Nullable
      column( :sensoroperationalmode, String ) #Nullable
      column( :beginposition, String )
      column( :procfacilityname, String ) #Nullable
      column( :orbitnumber, Int64 )
      column( :mode, String ) #Nullable
      column( :platformname, String )
      column( :producttype, String )
      column( :orbitdirection, String )
      column( :format, String )
      column( :passdirection, String ) #Nullable
      column( :procfacilityorg, String ) #Nullable
      column( :relativeorbitnumber, Int64 )
      column( :pduduration, Int64 ) #Nullable
      column( :platformidentifier, String )
      column( :gmlfootprint, String )
      column( :endposition, String )
      column( :creationdate, String ) #Nullable
      column( :onlinequalitycheck, String )
      column( :processingname, String ) #Nullable
      column( :identifier, String )
      column( :passnumber, Int64 ) #Nullable
      column( :pdualongtrackcoord, Int64 ) #Nullable
      column( :instrumentshortname, String )
      column( :processinglevel, String ) #Nullable
      column( :cloudcoverpercentage, Float64 ) #Nullable
      column( :lrmpercentage, Int64 ) #Nullable
      column( :openseapercentage, Int64 ) #Nullable
      column( :landpercentage, Int64 ) #Nullable
      column( :sarpercentage, Int64 ) #Nullable
      column( :closedseapercentage, Int64 ) #Nullable
      column( :continentalicepercentage, Int64 ) #Nullable
      column( :nbfire, Int64 ) #Nullable
      column( :lastpassdirection, String ) #Nullable
      column( :lastrelpassnumber, Int64 ) #Nullable
      column( :lastrelorbitdirection, String ) #Nullable
      column( :lastorbitdirection, String ) #Nullable
      column( :lastrelativeorbitnumber, Int64 ) #Nullable
      column( :lastpassnumber, Int64 ) #Nullable
      column( :lastrelpassdirection, String ) #Nullable
      column( :lastorbitnumber, Int64 ) #Nullable
      column( :pdutileid, String ) #Nullable
      column( :productclass, String ) #Nullable
      column( :polarisationmode, String ) #Nullable
      column( :status, String ) #Nullable
      column( :slicenumber, Int64 ) #Nullable
      column( :missiondatatakeid, Int64 ) #Nullable
      column( :swathidentifier, String ) #Nullable
      column( :acquisitiontype, String ) #Nullable
      column( :productconsolidation, String ) #Nullable
      column( :unclassifiedpercentage, Float64 ) #Nullable
      column( :s2datatakeid, String ) #Nullable
      column( :granuleidentifier, String ) #Nullable
      column( :level1cpdiidentifier, String ) #Nullable
      column( :vegetationpercentage, Float64 ) #Nullable
      column( :notvegetatedpercentage, Float64 ) #Nullable
      column( :processingbaseline, Float64 ) #Nullable
      column( :mediumprobacloudspercentage, Float64 ) #Nullable
      column( :datastripidentifier, String ) #Nullable
      column( :generationdate, String ) #Nullable
      column( :illuminationzenithangle, Float64 ) #Nullable
      column( :platformserialidentifier, String ) #Nullable
      column( :waterpercentage, Float64 ) #Nullable
      column( :highprobacloudspercentage, Float64 ) #Nullable
      column( :illuminationazimuthangle, Float64 ) #Nullable
      column( :snowicepercentage, Float64 ) #Nullable
      column( :hv_order_tileid, String ) #Nullable
      column( :tileid, String ) #Nullable
      column( :datatakesensingstart, String ) #Nullable
      column( :leapsecond, Int64 ) #Nullable
      column( :leapSecondOccurrence, String ) #Nullable
=#


#=
 uuid::String = ""
  relpassnumber::Int64 = 0 #Nullable
  relpassdirection::String = "" #Nullable
  filename::String = ""
  instrumentname::String = ""
  timeliness::String = "" #Nullable
  footprint::String = ""
  size::String = ""
  productlevel::String = "" #Nullable
  relorbitdir::String = "" #Nullable
  ingestiondate::String = ""
  ecmwf::String = "" #Nullable
  sensoroperationalmode::String = "" #Nullable
  beginposition::String = ""
  procfacilityname::String = "" #Nullable
  orbitnumber::Int64 = 0
  mode::String = "" #Nullable
  platformname::String = ""
  producttype::String = ""
  orbitdirection::String = ""
  format::String = ""
  passdirection::String = "" #Nullable
  procfacilityorg::String = "" #Nullable
  relativeorbitnumber::Int64 = 0
  pduduration::Int64 = 0 #Nullable
  platformidentifier::String = ""
  gmlfootprint::String = ""
  endposition::String = ""
  creationdate::String = "" #Nullable
  onlinequalitycheck::String = ""
  processingname::String = "" #Nullable
  identifier::String = ""
  passnumber::Int64 = 0 #Nullable
  pdualongtrackcoord::Int64 = 0 #Nullable
  instrumentshortname::String = ""
  processinglevel::String = "" #Nullable
  cloudcoverpercentage::Float64 = 0.0 #Nullable
  lrmpercentage::Int64 = 0 #Nullable
  openseapercentage::Int64 = 0 #Nullable
  landpercentage::Int64 = 0 #Nullable
  sarpercentage::Int64 = 0 #Nullable
  closedseapercentage::Int64 = 0 #Nullable
  continentalicepercentage::Int64 = 0 #Nullable
  nbfire::Int64 = 0 #Nullable
  lastpassdirection::String = "" #Nullable
  lastrelpassnumber::Int64 = 0 #Nullable
  lastrelorbitdirection::String = "" #Nullable
  lastorbitdirection::String = "" #Nullable
  lastrelativeorbitnumber::Int64 = 0 #Nullable
  lastpassnumber::Int64 = 0 #Nullable
  lastrelpassdirection::String = "" #Nullable
  lastorbitnumber::Int64 = 0 #Nullable
  pdutileid::String = "" #Nullable
  productclass::String = "" #Nullable
  polarisationmode::String = "" #Nullable
  status::String = "" #Nullable
  slicenumber::Int64 = 0 #Nullable
  missiondatatakeid::Int64 = 0 #Nullable
  swathidentifier::String = "" #Nullable
  acquisitiontype::String = "" #Nullable
  productconsolidation::String = "" #Nullable
  unclassifiedpercentage::Float64 = 0.0 #Nullable
  s2datatakeid::String = "" #Nullable
  granuleidentifier::String = "" #Nullable
  level1cpdiidentifier::String = "" #Nullable
  vegetationpercentage::Float64 = 0.0 #Nullable
  notvegetatedpercentage::Float64 = 0.0 #Nullable
  processingbaseline::Float64 = 0.0 #Nullable
  mediumprobacloudspercentage::Float64 = 0.0 #Nullable
  datastripidentifier::String = "" #Nullable
  generationdate::String = "" #Nullable
  illuminationzenithangle::Float64 = 0.0 #Nullable
  platformserialidentifier::String = "" #Nullable
  waterpercentage::Float64 = 0.0 #Nullable
  highprobacloudspercentage::Float64 = 0.0 #Nullable
  illuminationazimuthangle::Float64 = 0.0 #Nullable
  snowicepercentage::Float64 = 0.0 #Nullable
  hv_order_tileid::String = "" #Nullable
  tileid::String = "" #Nullable
  datatakesensingstart::String = "" #Nullable
  leapsecond::Int64 = 0 #Nullable
  leapSecondOccurrence::String = "" #Nullable

=#
