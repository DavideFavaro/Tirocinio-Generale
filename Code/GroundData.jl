module GroundData


using HTTP

using CombinedParsers
using CombinedParsers.Regexp

using ArchGDAL
using DataFrames
using Dates
using CSV

using Revise



@syntax informations = Sequence(
                            :begin => Sequence(
                                        re"<INIZIO>",
                                        Numeric(Int),
                                        re"</INIZIO>",
                                      ),
                            :end => Sequence(
                                        re"<FINE>",
                                        Numeric(Int),
                                        re"</FINE>",
                                    ),
                            :proj => Sequence(
                                        "<PROJECTION>",
                                        re"[^<:>]+",
                                        ":",
                                        Numeric(Int),
                                        "</PROJECTION>"
                                     )
                       )

@syntax stations = Repeat(
                        Sequence(
                            re"<[^<>]+>",
                            Either(
                                Numeric(Int),
                                Numeric(Float64),
                                Sequence(
                                    "<![CDATA[",
                                    re"[^\[\]]+",
                                    "]]>"
                                ),
                                re"[^<>]+"
                            ),
                            re"</[^<>]+>"
                        )
                    )


function getStationsInfo()

#Get the information on the available stations
    # Get the page containing informations on all the available stations
    page = HTTP.get( "https://www.arpa.veneto.it/bollettini/meteo60gg/stazioni.xml" )
    str_page = String(page.body)
    # Keep only the useful portion of the page
    useful = split( str_page, "</PERIODO>" )[2][1:end-15]

    # Split the downloaded string in substrings one for each station
    arr = split( useful, r"</?STAZIONE>", keepempty=false )

    info = informations(arr[1])
    stats = stations.(arr[2:end])

    dict = [ Dict( Symbol(String(attribute[1][2])) => attribute[2] for attribute in station ) for station in stats ]

    

    df = DataFrame(dict)

    return df
end


df = getStationsInfo()


end # moduled