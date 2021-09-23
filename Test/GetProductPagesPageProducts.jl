function createDict()
    io = getProductsPages( authenticate("davidefavaro","Tirocinio"), 200 )
    dict = getPageProducts( io )
    return dict
end