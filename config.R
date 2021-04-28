message("Setting system color palettes")
## pal1 = wes_palette("Royal1") #length 4
## pal2 = wes_palette("Zissou1") #length 5
## pal3 = wes_palette("IsleofDogs1") #length 6
## pal4 = wes_palette("Darjeeling1") #length 5
## pal5 = wes_palette("Moonrise1") #length 4
## pal6 = wes_palette("BottleRocket2") #length 5
## pal7 = wes_palette("Moonrise2") #length 4
## all.pal = c(pal1,pal2,pal3,pal4,pal5,pal6,pal7)

## xtYao's colors, not all will be used
etypes = c("del", "dup", "inv", "tra", "tic", "pyrgo", "bfb", "dm", "cpxdm",
           "chromoplexy", "chromothripsis", "invdup", "tyfonas", "rigma",
           "unclassified")

trubetskoy =
    setNames(c(
        "#e6194B", "#3cb44b", "#ffe119", "#4363d8", "#f58231", "#911eb4", "#42d4f4", "#f032e6", "#bfef45", "#fabebe", "#469990",
        "#e6beff", "#9A6324", "#fffac8", "#800000", "#aaffc3", "#808000", "#ffd8b1", "#000075", "#a9a9a9", "#ffffff", "#000000"),
        c("red", "green", "yellow", "blue", "orange", "purple", "cyan", "magenta", "lime", "pink", "teal",
          "lavender", "brown", "beige", "maroon", "mint", "olive", "apricot", "navy", "grey", "white", "black"))

event.cols = data.table(
    event = c(
        "del", "rigma", "chromothripsis", ## blue
        "qrp", "tic", ## green
        "tra", "chromoplexy", ## purple
        "inv", "invdup", ## yellow, because brown and 
        "dup", "pyrgo", ## orange
        "dm", "bfb", "tyfonas", ## red
        "other"
    ),
    color = c("blue", "navy", "steel",
              "teal", "green",
              "olive", "lime",
              "magenta", "purple",
              "red", "maroon",
              "yellow", "orange", "brown",
              "grey50")
)
