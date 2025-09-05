################################################################################
## Script Name:        Format Seabird Orders
## Purpose:            Regex functions for formatting seabird names (using birdlife)
##                      taxonomy
## Author:             James Baxter
## Date Created:       2025-09-05
################################################################################

############################### SYSTEM OPTIONS #################################
options(
  scipen = 6,     # Avoid scientific notation
  digits = 7      # Set precision for numerical display
)
memory.limit(30000000)

############################### DEPENDENCIES ###################################
# Load required libraries
library(tidyverse)
library(magrittr)


################################## Functions ###################################
# Read and inspect data
FormatCharadriiformes <- function(x) {
  charadriiformes <- c(
    "Chionis albus" = "Snowy Sheathbill|Chionis albus",
    "Chionis minor" = "(B|b)lack-faced Sheathbill|Chionis minor",
    "Anous stolidus" = "(B|b)rown Noddy|Anous stolidus",
    "Anous tenuirostris" = "Lesser Noddy|Anous tenuirostris",
    "Anous minutus" = "(B|b)lack Noddy|Anous minutus",
    "Anous ceruleus" = "(B|b)lue Noddy|Anous ceruleus",
    "Anous albivittus" = "Grey Noddy|Anous albivittus",
    "Gygis alba" = "Atlantic White[ -_](T|t)ern|Gygis alba",
    "Gygis candida" = "Common White[ -_](T|t)ern|Gygis candida",
    "Gygis microrhyncha" = "Little White[ -_](T|t)ern|Gygis microrhyncha",
    "Rynchops niger" = "(B|b)lack[ _-](S|s)kimmer|Rynchops niger",
    "Rynchops flavirostris" = "African Skimmer|Rynchops flavirostris",
    "Rynchops albicollis" = "Indian Skimmer|Rynchops albicollis",
    "Saundersilarus saundersi" = "Saunders('{0,1}s){0,1}[ -_](G|g)ull|Saundersilarus saundersi",
    "Hydrocoloeus minutus" = "(L|l)ittle[ -_](G|g)ull|Hydrocoloeus minutus",
    "Rhodostethia rosea" = "Ross('{0,1}s){0,1}[ -_](G|g)ull|Rhodostethia rosea",
    "Creagrus furcatus" = "Swallow-tailed[ -_](G|g)ull|Creagrus furcatus",
    "Xema sabini" = "Sabine('{0,1}s){0,1}[ -_](G|g)ull|Xema sabini",
    "Pagophila eburnea" = "Ivory[ -_](G|g)ull|Pagophila eburnea",
    "Rissa brevirostris" = "Red[ -_]legged[ -_](K|k)ittiwake|Rissa brevirostris",
    "Rissa tridactyla" = "(B|b)lack[ -_]legged[ -_](K|k)ittiwake|Rissa tridactyla",
    "Larus philadelphia" = "(B|b)onaparte('{0,1}s){0,1}[ -_](G|g)ull|(Larus|Chroicocephalus) philadelphia",
    "Larus genei" = "(S|s)lender[ -_]billed[ -_](G|g)ull|Larus genei",
    "Larus brunnicephalus" = "(B|b)rown-headed[ -_](G|g)ull|Larus brunnicephalus",
    "Larus ridibundus" = "(B|b)lac(k|h)[ -_](h|H)eadea{0,1}d[ -_](sea){0,1}(G|g)ull|(Larus|Chroicocephalus)[ -_]ridibundus|(B|b)l_(H|h)_gull",
    "Larus serranus" = "Andean[ -_](G|g)ull|Larus serranus",
    "Larus maculipennis" = "(B|b)rown-hooded[ -_](G|g)ull|Larus maculipennis",
    "Larus hartlaubii" = "Hartlaub('{0,1}s){0,1}[ -_](G|g)ull|Larus hartlaubii",
    "Larus cirrocephalus" = "Gr(e|a)y[ -_](h|H)eaded[ -_](G|g)ull|Larus[ -_]{0,1}cirrocephalus",
    "Larus novaehollandiae" = "Silver[ -_](G|g)ull|Larus novaehollandiae",
    "Larus bulleri" = "(B|b)lack-billed[ -_](G|g)ull|Larus bulleri",
    "Larus modestus" = "(G|g)r(e|a)y[ -_](G|g)ull|Larus modestus",
    "Larus scoresbii" = "Dolphin[ -_](G|g)ull|Larus scoresbii",
    "Larus pipixcan" = "Franklin('{0,1}s){0,1}[ -_](G|g)ull|Larus pipixcan",
    "Larus atricilla" = "Laughing[ -_](G|g)ull|Larus atricilla",
    "Larus fuliginosus" = "Lava[ -_](G|g)ull|Larus fuliginosus",
    "Larus ichthyaetus" = "Pallas('{0,1}s){0,1}[ -_](G|g)ull|Larus ichthyaetus",
    "Larus relictus" = "Relict[ -_](G|g)ull|Larus relictus",
    "Larus melanocephalus" = "(M|m)editerranean[ -_](G|g)ull|(Larus|Ichthyaetus) melanocephalus",
    "Larus hemprichii" = "Sooty[ -_](G|g)ull|Larus hemprichii",
    "Larus leucophthalmus" = "White-eyed[ -_](G|g)ull|Larus leucophthalmus",
    "Larus audouinii" = "Audouin('{0,1}s){0,1}[ -_](G|g)ull|Larus audouinii",
    "Larus heermanni" = "Heermann('{0,1}s){0,1}[ -_](G|g)ull|Larus heermanni",
    "Larus pacificus" = "Pacific[ -_](G|g)ull|Larus pacificus",
    "Larus crassirostris" = "(B|b)lack-tailed[ -_](G|g)ull|Larus crassirostris",
    "Larus belcheri" = "((B|b)elcher('{0,1}s){0,1}|(b|B)and[ -_](t|T)ailed)[ -_](G|g)ull|Larus belcheri",
    "Larus atlanticus" = "Olrog('{0,1}s){0,1}[ -_](G|g)ull|Larus atlanticus",
    "Larus delawarensis" = "(R|r)ing[ -_](B|b)illed[ -_](G|g)ull|Larus delawarensis",
    "Larus canus" = "Mew[ -_](G|g)ull|Larus[ -_]canus|(s|S)hort[ -_](b|B)illed[ -_](G|g)ull|(C|c)ommon[ -_]{0,1}(G|g)ull|(c|C)[ -_](G|g)ull",
    "Larus livens" = "Yellow-footed[ -_](G|g)ull|Larus[ -_]livens",
    "Larus occidentalis" = "(W|w)estern (G|g)ull|Larus occidentalis",
    "Larus californicus" = "California[ -_](G|g)ull|Larus californicus",
    "Larus dominicanus" = "(K|k)elp[ _-](G|g)ull|Larus dominicanus",
    "Larus fuscus" = "(L|l)esser[ _-](B|b)lack[- ](B|b)acked[ _-](G|g)ull|Larus fuscus|(L|l)[ -_](B|b)l[ -_](B|b)a[ -_]{0,1}(G|g)ull",
    "Larus argentatus" = "((E|e)uropean[ _-]{0,1}){0,1}(H|h)erring{0,1}[ _-]{0,1}(G|g)ull|Larus[ -_]{0,1}argentatus",
    "Larus armenicus" = "Armenian[ -_](G|g)ull|Larus armenicus",
    "Larus michahellis" = "(Y|y)ellow[ _-]{0,1}(L|l)eggen{0,1}d[ _-]{0,1}(G|g)ull|Larus[ _-]michah{0,1}ellis",
    "Larus cachinnans" = "Caspian[ -_](G|g)ull|Larus cachinnans",
    "Larus smithsonianus" = "Arctic Herring[ -_](G|g)ull|Larus smithsonianus",
    "Larus glaucoides" = "(Iceland|(t|T)hayers)[ -_](G|g)ull|Larus glaucoides",
    "Larus schistisagus" = "(s|S)laty[ -_]{0,1}(B|b)acked[ _-](G|g)ull|Larus schistisagus",
    "Larus glaucescens" = "Glaucous-winged[ -_](G|g)ull|Larus glaucescens",
    "Larus hyperboreus" = "(G|g)laucous[ _-](G|g)ull|Larus hyperboreus",
    "Larus marinus" = "(G|g)reat[ _-](B|b)lack[ _-](B|b)ac{0,1}ked[ _-](G|g)ull|Larus[ -_]marinus|(G|g)r[ _-](B|b)k[ _-](B|b)d[ _-](G|g)ull",
    "Onychoprion aleuticus" = "Aleutian[ -_](T|t)ern|Onychoprion aleuticus",
    "Onychoprion fuscatus" = "Sooty[ -_](T|t)ern|Onychoprion fuscatus",
    "Onychoprion anaethetus" = "(B|b)ridled[ -_](T|t)ern|Onychoprion anaethetus",
    "Onychoprion lunatus" = "Grey-backed[ -_](T|t)ern|Onychoprion lunatus",
    "Sternula albifrons" = "(L|l)ittle[ -_](T|t)ern|Sternula albifrons",
    "Sternula saundersi" = "Saunders('{0,1}s){0,1}[ -_](T|t)ern|Sternula saundersi",
    "Sternula antillarum" = "Least[ -_](T|t)ern|Sternula antillarum",
    "Sternula superciliaris" = "Yellow-billed[ -_](T|t)ern|Sternula superciliaris",
    "Sternula lorata" = "Peruvian[ -_](T|t)ern|Sternula lorata",
    "Sternula nereis" = "Fairy[ -_](T|t)ern|Sternula nereis",
    "Sternula balaenarum" = "Damara[ -_](T|t)ern|Sternula balaenarum",
    "Phaetusa simplex" = "Large-billed[ -_](T|t)ern|Phaetusa simplex",
    "Gelochelidon nilotica" = "(G|g)ull-billed[ -_](T|t)ern|Gelochelidon nilotica",
    "Hydroprogne caspia" = "(c|C)aspian (T|t)ern|Hydroprogne caspia",
    "Larosterna inca" = "Inca[ -_](T|t)ern|Larosterna inca",
    "Chlidonias albostriatus" = "(B|b)lack-fronted[ -_](T|t)ern|Chlidonias albostriatus",
    "Chlidonias hybrida" = "(W|w)hiskered[ -_](T|t)ern|Chlidonias hybrida",
    "Chlidonias leucopterus" = "White-winged[ -_](T|t)ern|Chlidonias leucopterus",
    "Chlidonias niger" = "(B|b)lack[ -_](T|t)ern|Chlidonias niger",
    "Sterna aurantia" = "River[ -_](T|t)ern|Sterna aurantia",
    "Sterna dougallii" = "Roseate[ -_](T|t)ern|Sterna dougallii",
    "Sterna striata" = "White-fronted[ -_](T|t)ern|Sterna striata",
    "Sterna sumatrana" = "(B|b)lack-naped[ -_](T|t)ern|Sterna sumatrana",
    "Sterna hirundinacea" = "(s|S)outh[ -_]{0,1}(A|a)merican[ -_](T|t)ern|Sterna hirundinacea",
    "Sterna hirundo" = "(c|C)ommon[ -_]{0,1}(t|T)ern|Sterna[ _-]hirundo",
    "Sterna repressa" = "White-cheeked[ -_](T|t)ern|Sterna repressa",
    "Sterna paradisaea" = "(A|a)rc{0,1}tic[ -_](T|t)ern|Sterna paradisaea",
    "Sterna vittata" = "Antarctic[ -_](T|t)ern|Sterna vittata",
    "Sterna virgata" = "Kerguelen[ -_](T|t)ern|Sterna virgata",
    "Sterna forsteri" = "Forster('{0,1}s){0,1}[ -_](T|t)ern|Sterna forsteri",
    "Sterna trudeaui" = "Snowy-crowned[ -_](T|t)ern|Sterna trudeaui",
    "Sterna acuticauda" = "(B|b)lack-bellied[ -_](T|t)ern|Sterna acuticauda",
    "Thalasseus bengalensis" = "Lesser Crested[ -_](T|t)ern|Thalasseus bengalensis",
    "Thalasseus bernsteini" = "Chinese Crested[ -_](T|t)ern|Thalasseus bernsteini",
    "Thalasseus elegans" = "Elegant[ -_](T|t)ern|Thalasseus elegans", 
    "Thalasseus sandvicensis" = "(s|S)andwich[ -_](T|t)ern|(Thalasseus|Sterna) sandvicensis",
    "Thalasseus maximus" = "(R|r)oyal[ -_](T|t)ern|Thalasseus maximus",
    "Thalasseus bergii" = "(g|G)reater[ -_](c|C)rested[ -_](T|t)ern|Thalasseus bergii|(S|s)wift[ -_]{0,1}(T|t)ern",
    "Stercorarius longicaudus" = "(L|l)ong[ _-](T|t)ailed[ _-]((j|J)aeger|(s|S)kua)|Stercorarius longicaudus",
    "Stercorarius parasiticus" = "Arctic Jaeger|Stercorarius parasiticus",
    "Stercorarius pomarinus" = "Pomarine Jaeger|Stercorarius pomarinus",
    "Catharacta skua" = "(G|g)reat[ -_](S|s)kua|(Catharacta|Stercorarius) skua",
    "Catharacta maccormicki" = "(S|s)outh[ -_](P|p)olar[ -_](S|s)kua|(Catharacta|Stercorarius) maccormicki",
    "Catharacta antarctica" = "(B|b)rown[ -_](S|s)kua|(Catharacta|Stercorarius) antarctica",
    "Catharacta chilensis" = "Chilean Skua|(Catharacta|Stercorarius) chilensis",
    "Cerorhinca monocerata" = "Rhinoceros Auklet|(Catharacta|Stercorarius) monocerata",
    "Fratercula cirrhata" = "Tufted Puffin|Fratercula cirrhata",
    "Fratercula arctica" = "(A|a)tlantic[ -_](P|p)uffin|Fratercula arctica",
    "Fratercula corniculata" = "Horned Puffin|Fratercula corniculata",
    "Ptychoramphus aleuticus" = "Cassin('{0,1}s){0,1} Auklet|Ptychoramphus aleuticus",
    "Aethia psittacula" = "Parakeet Auklet|Aethia psittacula",
    "Aethia pusilla" = "Least Auklet|Aethia pusilla",
    "Aethia pygmaea" = "Whiskered Auklet|Aethia pygmaea",
    "Aethia cristatella" = "Crested Auklet|Aethia cristatella",
    "Brachyramphus perdix" = "Long-billed Murrelet|Brachyramphus perdix",
    "Brachyramphus marmoratus" = "Marbled Murrelet|Brachyramphus marmoratus",
    "Brachyramphus brevirostris" = "Kittlitz('{0,1}s){0,1} Murrelet|Brachyramphus brevirostris",
    "Cepphus grylle" = "(B|b)lack Guillemot|Cepphus grylle",
    "Cepphus columba" = "Pigeon Guillemot|Cepphus columba",
    "Cepphus carbo" = "Spectacled Guillemot|Cepphus carbo",
    "Synthliboramphus antiquus" = "Ancient Murrelet|Synthliboramphus antiquus",
    "Synthliboramphus wumizusume" = "Japanese Murrelet|Synthliboramphus wumizusume",
    "Synthliboramphus scrippsi" = "Scripps('{0,1}s){0,1} Murrelet|Synthliboramphus scrippsi",
    "Synthliboramphus hypoleucus" = "Guadalupe Murrelet|Synthliboramphus hypoleucus",
    "Synthliboramphus craveri" = "Craveri('{0,1}s){0,1} Murrelet|Synthliboramphus craveri",
    "Alca torda" = "Razorbill|Alca torda",
    "Pinguinus impennis" = "Great Auk|Pinguinus impennis",
    "Alle alle" = "Little Auk|Alle alle",
    "Uria lomvia" = "(T|t)hick[ -_]billed[ -_](M|m)urre|Uria lomvia",
    "Uria aalge" = "(C|c)ommon[ -_]((M|m)urre|(g|G)uillemot)|Uria aalge",
    "Laridae sp." = '^(sea){0,1}(G|g)ull\\d{0,3}$|^(l|L)arus$|^(t|T)ern$|^(b|b)lack[ -_]backed[ -_]gull',
    'Alcidae sp.' = '^(g|G)uillemot$'
  )
  
  z <- NA_character_
  
  for (i in 1:length(charadriiformes)) {
    if (any(grepl(charadriiformes[[i]], x))) {
      z <- names(charadriiformes)[[i]]
    }
  }
  
  return(z)
}


FormatProcellariiformes <- function(x) {
  procellariiformes <- c(
    "Oceanites oceanicus" = "Wilson('{0,1}s){0,1} Storm-petrel|Oceanites oceanicus",
    "Oceanites gracilis" = "White-vented Storm-petrel|Oceanites gracilis",
    "Oceanites pincoyae" = "Pincoya Storm-petrel|Oceanites pincoyae",
    "Garrodia nereis" = "Grey-backed Storm-petrel|Garrodia nereis",
    "Pelagodroma marina" = "White-faced Storm-petrel|Pelagodroma marina",
    "Fregetta grallaria" = "White-bellied Storm-petrel|Fregetta grallaria",
    "Fregetta tropica" = "(B|b)lack-bellied Storm-petrel|Fregetta tropica",
    "Fregetta maoriana" = "New Zealand Storm-petrel|Fregetta maoriana",
    "Fregetta lineata" = "New Caledonian Storm-petrel|Fregetta lineata",
    "Nesofregetta fuliginosa" = "Polynesian Storm-petrel|Nesofregetta fuliginosa",
    "Hydrobates pelagicus" = "European Storm-petrel|Hydrobates pelagicus",
    "Hydrobates jabejabe" = "Cape Verde Storm-petrel|Hydrobates jabejabe",
    "Hydrobates castro" = "(B|b)and-rumped Storm-petrel|Hydrobates castro",
    "Hydrobates monteiroi" = "Monteiro('{0,1}s){0,1} Storm-petrel|Hydrobates monteiroi",
    "Hydrobates matsudairae" = "Matsudaira('{0,1}s){0,1} Storm-petrel|Hydrobates matsudairae",
    "Hydrobates melania" = "(B|b)lack Storm-petrel|Hydrobates melania",
    "Hydrobates homochroa" = "Ashy Storm-petrel|Hydrobates homochroa",
    "Hydrobates microsoma" = "Least Storm-petrel|Hydrobates microsoma",
    "Hydrobates tethys" = "Wedge-rumped Storm-petrel|Hydrobates tethys",
    "Hydrobates socorroensis" = "Townsend('{0,1}s){0,1} Storm-petrel|Hydrobates socorroensis",
    "Hydrobates cheimomnestes" = "Ainley('{0,1}s){0,1} Storm-petrel|Hydrobates cheimomnestes",
    "Hydrobates leucorhous" = "Leach('{0,1}s){0,1} Storm-petrel|Hydrobates leucorhous",
    "Hydrobates monorhis" = "Swinhoe('{0,1}s){0,1} Storm-petrel|Hydrobates monorhis",
    "Hydrobates macrodactylus" = "Guadalupe Storm-petrel|Hydrobates macrodactylus",
    "Hydrobates tristrami" = "Tristram('{0,1}s){0,1} Storm-petrel|Hydrobates tristrami",
    "Hydrobates markhami" = "Markham('{0,1}s){0,1} Storm-petrel|Hydrobates markhami",
    "Hydrobates furcatus" = "Fork-tailed Storm-petrel|Hydrobates furcatus",
    "Hydrobates hornbyi" = "Ringed Storm-petrel|Hydrobates hornbyi",
    "Diomedea sanfordi" = "Northern Royal Albatross|Diomedea sanfordi",
    "Diomedea epomophora" = "Southern Royal Albatross|Diomedea epomophora",
    "Diomedea exulans" = "Snowy Albatross|Diomedea exulans",
    "Diomedea antipodensis" = "Antipodean Albatross|Diomedea antipodensis",
    "Diomedea amsterdamensis" = "Amsterdam Albatross|Diomedea amsterdamensis",
    "Diomedea dabbenena" = "Tristan Albatross|Diomedea dabbenena",
    "Phoebetria fusca" = "Sooty Albatross|Phoebetria fusca",
    "Phoebetria palpebrata" = "Light-mantled Albatross|Phoebetria palpebrata",
    "Phoebastria irrorata" = "Waved Albatross|Phoebastria irrorata",
    "Phoebastria nigripes" = "(B|b)lack-footed Albatross|Phoebastria nigripes",
    "Phoebastria immutabilis" = "Laysan Albatross|Phoebastria immutabilis",
    "Phoebastria albatrus" = "Short-tailed Albatross|Phoebastria albatrus",
    "Thalassarche chlororhynchos" = "Atlantic Yellow-nosed Albatross|Thalassarche chlororhynchos",
    "Thalassarche carteri" = "Indian Yellow-nosed Albatross|Thalassarche carteri",
    "Thalassarche chrysostoma" = "Grey-headed Albatross|Thalassarche chrysostoma",
    "Thalassarche melanophris" = "(B|b)lack-browed Albatross|Thalassarche melanophris",
    "Thalassarche impavida" = "Campbell Albatross|Thalassarche impavida",
    "Thalassarche bulleri" = "(B|b)uller('{0,1}s){0,1} Albatross|Thalassarche bulleri",
    "Thalassarche cauta" = "Shy Albatross|Thalassarche cauta",
    "Thalassarche steadi" = "White-capped Albatross|Thalassarche steadi",
    "Thalassarche eremita" = "Chatham Albatross|Thalassarche eremita",
    "Thalassarche salvini" = "Salvin('{0,1}s){0,1} Albatross|Thalassarche salvini",
    "Macronectes halli" = "(N|n)orthern[ -_](G|g)iant[ -_](P|p)etrel|Macronectes halli",
    "Macronectes giganteus" = "Southern Giant[ -_]{0,1}(P|p)etrel|Macronectes giganteus",
    "Fulmarus glacialis" = "(N|n)orthern[ -_](F|f)ulmar|Fulmarus glacialis",
    "Fulmarus glacialoides" = "(S|s)outhern[ -_](F|f)ulmar|Fulmarus glacialoides",
    "Thalassoica antarctica" = "Antarctic[ -_]{0,1}(P|p)etrel|Thalassoica antarctica",
    "Daption capense" = "Cape[ -_]{0,1}(P|p)etrel|Daption capense",
    "Pagodroma nivea" = "Snow[ -_]{0,1}(P|p)etrel|Pagodroma nivea",
    "Halobaena caerulea" = "(B|b)lue[ -_]{0,1}(P|p)etrel|Halobaena caerulea",
    "Pachyptila vittata" = "(B|b)road-billed Prion|Pachyptila vittata",
    "Pachyptila salvini" = "Salvin('{0,1}s){0,1} Prion|Pachyptila salvini",
    "Pachyptila macgillivrayi" = "MacGillivray('{0,1}s){0,1} Prion|Pachyptila macgillivrayi",
    "Pachyptila desolata" = "Antarctic Prion|Pachyptila desolata",
    "Pachyptila belcheri" = "Slender-billed Prion|Pachyptila belcheri",
    "Pachyptila turtur" = "Fairy Prion|Pachyptila turtur",
    "Pachyptila crassirostris" = "Fulmar Prion|Pachyptila crassirostris",
    "Aphrodroma brevirostris" = "Kerguelen[ -_]{0,1}(P|p)etrel|Aphrodroma brevirostris",
    "Pterodroma rupinarum" = "Large Saint Helena[ -_]{0,1}(P|p)etrel|Pterodroma rupinarum",
    "Pterodroma leucoptera" = "White-winged[ -_]{0,1}(P|p)etrel|Pterodroma leucoptera",
    "Pterodroma brevipes" = "Collared[ -_]{0,1}(P|p)etrel|Pterodroma brevipes",
    "Pterodroma defilippiana" = "Masatierra[ -_]{0,1}(P|p)etrel|Pterodroma defilippiana",
    "Pterodroma longirostris" = "Stejneger('{0,1}s){0,1}[ -_]{0,1}(P|p)etrel|Pterodroma longirostris",
    "Pterodroma cookii" = "Cook('{0,1}s){0,1}[ -_]{0,1}(P|p)etrel|Pterodroma cookii",
    "Pterodroma pycrofti" = "Pycroft('{0,1}s){0,1}[ -_]{0,1}(P|p)etrel|Pterodroma pycrofti",
    "Pterodroma hypoleuca" = "(B|b)onin[ -_]{0,1}(P|p)etrel|Pterodroma hypoleuca",
    "Pterodroma nigripennis" = "(B|b)lack-winged[ -_]{0,1}(P|p)etrel|Pterodroma nigripennis",
    "Pterodroma axillaris" = "Chatham Islands[ -_]{0,1}(P|p)etrel|Pterodroma axillaris",
    "Pterodroma ultima" = "Murphy('{0,1}s){0,1}[ -_]{0,1}(P|p)etrel|Pterodroma ultima",
    "Pterodroma solandri" = "Providence[ -_]{0,1}(P|p)etrel|Pterodroma solandri",
    "Pterodroma neglecta" = "Kermadec[ -_]{0,1}(P|p)etrel|Pterodroma neglecta",
    "Pterodroma arminjoniana" = "Trindade[ -_]{0,1}(P|p)etrel|Pterodroma arminjoniana",
    "Pterodroma heraldica" = "Herald[ -_]{0,1}(P|p)etrel|Pterodroma heraldica",
    "Pterodroma atrata" = "Henderson[ -_]{0,1}(P|p)etrel|Pterodroma atrata",
    "Pterodroma alba" = "Phoenix[ -_]{0,1}(P|p)etrel|Pterodroma alba",
    "Pterodroma baraui" = "(B|b)arau('{0,1}s){0,1}[ -_]{0,1}(P|p)etrel|Pterodroma baraui",
    "Pterodroma inexpectata" = "Mottled[ -_]{0,1}(P|p)etrel|Pterodroma inexpectata",
    "Pterodroma sandwichensis" = "Hawaiian[ -_]{0,1}(P|p)etrel|Pterodroma sandwichensis",
    "Pterodroma phaeopygia" = "Galapagos[ -_]{0,1}(P|p)etrel|Pterodroma phaeopygia",
    "Pterodroma cervicalis" = "White-necked[ -_]{0,1}(P|p)etrel|Pterodroma cervicalis",
    "Pterodroma externa" = "Juan Fernandez[ -_]{0,1}(P|p)etrel|Pterodroma externa",
    "Pterodroma mollis" = "Soft-plumaged[ -_]{0,1}(P|p)etrel|Pterodroma mollis",
    "Pterodroma cahow" = "(B|b)ermuda[ -_]{0,1}(P|p)etrel|Pterodroma cahow",
    "Pterodroma hasitata" = "(B|b)lack-capped[ -_]{0,1}(P|p)etrel|Pterodroma hasitata",
    "Pterodroma caribbaea" = "Jamaican[ -_]{0,1}(P|p)etrel|Pterodroma caribbaea",
    "Pterodroma feae" = "Cape Verde[ -_]{0,1}(P|p)etrel|Pterodroma feae",
    "Pterodroma deserta" = "Desertas[ -_]{0,1}(P|p)etrel|Pterodroma deserta",
    "Pterodroma madeira" = "Zino('{0,1}s){0,1}[ -_]{0,1}(P|p)etrel|Pterodroma madeira",
    "Pterodroma magentae" = "Magenta[ -_]{0,1}(P|p)etrel|Pterodroma magentae",
    "Pterodroma incerta" = "Atlantic[ -_]{0,1}(P|p)etrel|Pterodroma incerta",
    "Pterodroma lessonii" = "White-headed[ -_]{0,1}(P|p)etrel|Pterodroma lessonii",
    "Pterodroma macroptera" = "Great-winged[ -_]{0,1}(P|p)etrel|Pterodroma macroptera",
    "Pterodroma gouldi" = "Grey-faced[ -_]{0,1}(P|p)etrel|Pterodroma gouldi",
    "Procellaria cinerea" = "Grey[ -_]{0,1}(P|p)etrel|Procellaria cinerea",
    "Procellaria aequinoctialis" = "White-chinned[ -_]{0,1}(P|p)etrel|Procellaria aequinoctialis",
    "Procellaria conspicillata" = "Spectacled[ -_]{0,1}(P|p)etrel|Procellaria conspicillata",
    "Procellaria westlandica" = "Westland[ -_]{0,1}(P|p)etrel|Procellaria westlandica",
    "Procellaria parkinsoni" = "(B|b)lack[ -_]{0,1}(P|p)etrel|Procellaria parkinsoni",
    "Ardenna pacifica" = "Wedge-tailed[ -_](S|s)hearwater|Ardenna pacifica",
    "Ardenna bulleri" = "(B|b)uller('{0,1}s){0,1}[ -_](S|s)hearwater|Ardenna bulleri",
    "Ardenna tenuirostris" = "(S|s)hort[ -_](t|T)ailed[ -_](S|s)hearwater|Ardenna tenuirostris",
    "Ardenna grisea" = "Sooty[ -_](S|s)hearwater|Ardenna grisea",
    "Ardenna gravis" = "(G|g)reat[ -_](S|s)hearwater|Ardenna gravis",
    "Ardenna carneipes" = "Flesh-footed[ -_](S|s)hearwater|Ardenna carneipes",
    "Ardenna creatopus" = "Pink-footed[ -_](S|s)hearwater|Ardenna creatopus",
    "Calonectris leucomelas" = "Streaked[ -_](S|s)hearwater|Calonectris leucomelas",
    "Calonectris diomedea" = "Scopoli('{0,1}s){0,1}[ -_](S|s)hearwater|Calonectris diomedea",
    "Calonectris borealis" = "Cory('{0,1}s){0,1}[ -_](S|s)hearwater|Calonectris borealis",
    "Calonectris edwardsii" = "Cape Verde[ -_](S|s)hearwater|Calonectris edwardsii",
    "Puffinus nativitatis" = "Christmas[ -_](S|s)hearwater|Puffinus nativitatis",
    "Puffinus subalaris" = "Galapagos[ -_](S|s)hearwater|Puffinus subalaris",
    "Puffinus gavia" = "Fluttering[ -_](S|s)hearwater|Puffinus gavia",
    "Puffinus huttoni" = "Hutton('{0,1}s){0,1}[ -_](S|s)hearwater|Puffinus huttoni",
    "Puffinus opisthomelas" = "(B|b)lack-vented[ -_](S|s)hearwater|Puffinus opisthomelas",
    "Puffinus bryani" = "(B|b)ryan('{0,1}s){0,1}[ -_](S|s)hearwater|Puffinus bryani",
    "Puffinus myrtae" = "Rapa[ -_](S|s)hearwater|Puffinus myrtae",
    "Puffinus newelli" = "Newell('{0,1}s){0,1}[ -_](S|s)hearwater|Puffinus newelli",
    "Puffinus auricularis" = "Townsend('{0,1}s){0,1}[ -_](S|s)hearwater|Puffinus auricularis",
    "Puffinus bailloni" = "Tropical[ -_](S|s)hearwater|Puffinus bailloni",
    "Puffinus persicus" = "Persian[ -_](S|s)hearwater|Puffinus persicus",
    "Puffinus bannermani" = "(B|b)annerman('{0,1}s){0,1}[ -_](S|s)hearwater|Puffinus bannermani",
    "Puffinus puffinus" = "(M|m)anx[ -_](S|s)hearwater|Puffinus puffinus",
    "Puffinus yelkouan" = "Yelkouan[ -_](S|s)hearwater|Puffinus yelkouan",
    "Puffinus mauretanicus" = "(B|b)alearic[ -_](S|s)hearwater|Puffinus mauretanicus",
    "Puffinus elegans" = "Subantarctic[ -_](S|s)hearwater|Puffinus elegans",
    "Puffinus assimilis" = "Little[ -_](S|s)hearwater|Puffinus assimilis",
    "Puffinus lherminieri" = "Audubon('{0,1}s){0,1}[ -_](S|s)hearwater|Puffinus lherminieri",
    "Puffinus heinrothi" = "Heinroth('{0,1}s){0,1}[ -_](S|s)hearwater|Puffinus heinrothi",
    "Pseudobulweria macgillivrayi" = "Fiji[ -_]{0,1}(P|p)etrel|Pseudobulweria macgillivrayi",
    "Pseudobulweria aterrima" = "Mascarene[ -_]{0,1}(P|p)etrel|Pseudobulweria aterrima",
    "Pseudobulweria becki" = "(B|b)eck('{0,1}s){0,1}[ -_]{0,1}(P|p)etrel|Pseudobulweria becki",
    "Pseudobulweria rostrata" = "Tahiti[ -_]{0,1}(P|p)etrel|Pseudobulweria rostrata",
    "(B|b)ulweria bulwerii" = "(B|b)ulwer('{0,1}s){0,1}[ -_]{0,1}(P|p)etrel|Bulweria bulwerii",
    "(B|b)ulweria fallax" = "Jouanin('{0,1}s){0,1}[ -_]{0,1}(P|p)etrel|Bulweria fallax",
    "(B|b)ulweria bifax" = "Small Saint Helena[ -_]{0,1}(P|p)etrel|Bulweria bifax",
    "Pelecanoides whenuahouensis" = "Whenua Hou Diving-petrel|Pelecanoides whenuahouensis",
    "Pelecanoides garnotii" = "Peruvian Diving-petrel|Pelecanoides garnotii",
    "Pelecanoides magellani" = "Magellanic Diving-petrel|Pelecanoides magellani",
    "Pelecanoides georgicus" = "South Georgia Diving-petrel|Pelecanoides georgicus",
    "Pelecanoides urinatrix" = "Common Diving-petrel|Pelecanoides urinatrix",
    'Procellariidae sp.' = '^(F|f)ulmar$'
  )
  
  z <- NA_character_
  
  for (i in 1:length(procellariiformes)) {
    if (any(grepl(procellariiformes[[i]], x))) {
      z <- names(procellariiformes)[[i]]
    }
  }
  
  return(z)
}


FormatSphenisciformes <- function(x) {
  sphenisciformes <- c(
    "Aptenodytes patagonicus" = "King Penguin|Aptenodytes patagonicus",
    "Aptenodytes forsteri" = "Emperor Penguin|Aptenodytes forsteri",
    "Pygoscelis papua" = "Gentoo Penguin|Pygoscelis papua",
    "Pygoscelis adeliae" = "Adelie Penguin|Pygoscelis adeliae",
    "Pygoscelis antarcticus" = "Chinstrap Penguin|Pygoscelis antarcticus",
    "Eudyptes schlegeli" = "Royal Penguin|Eudyptes schlegeli",
    "Eudyptes chrysolophus" = "Macaroni Penguin|Eudyptes chrysolophus",
    "Eudyptes moseleyi" = "Northern Rockhopper Penguin|Eudyptes moseleyi",
    "Eudyptes chrysocome" = "Southern Rockhopper Penguin|Eudyptes chrysocome",
    "Eudyptes sclateri" = "Erect-crested Penguin|Eudyptes sclateri",
    "Eudyptes pachyrhynchus" = "Fiordland Penguin|Eudyptes pachyrhynchus",
    "Eudyptes robustus" = "Snares Penguin|Eudyptes robustus",
    "Megadyptes antipodes" = "Yellow-eyed Penguin|Megadyptes antipodes",
    "Eudyptula minor" = "Little Penguin|Eudyptula minor",
    "Spheniscus demersus" = "(A|a)frican[ -_]{0,1}(P|p)enguin|Spheniscus demersus",
    "Spheniscus magellanicus" = "Magellanic Penguin|Spheniscus magellanicus",
    "Spheniscus humboldti" = "(H|h)umboldt[ -_]{0,1}(P|p)enguin|Spheniscus humboldti",
    "Spheniscus mendiculus" = "Galapagos Penguin|Spheniscus mendiculus",
    "Spheniscidae sp." = '^(p|P)enguin$'
  )
  
  z <- NA_character_
  
  for (i in 1:length(sphenisciformes)) {
    if (any(grepl(sphenisciformes[[i]], x))) {
      z <- names(sphenisciformes)[[i]]
    }
  }
  
  return(z)
}


FormatSuliformes <- function(x) {
  suliformes <- c(
    "Fregata ariel" = "Lesser Frigatebird|Fregata ariel",
    "Fregata minor" = "Great Frigatebird|Fregata minor",
    "Fregata andrewsi" = "Christmas Island Frigatebird|Fregata andrewsi",
    "Fregata magnificens" = "Magnificent Frigatebird|Fregata magnificens",
    "Fregata aquila" = "Ascension Frigatebird|Fregata aquila",
    "Papasula abbotti" = "Abbott('{0,1}s){0,1} Booby|Papasula abbotti",
    "Morus bassanus" = "(n|N)or{0,1}thern[ -_]{0,1}(g|G)annet|Morus bassanus",
    "Morus capensis" = "(C|c)ape (G|g)annet|Morus capensis",
    "Morus serrator" = "Australasian Gannet|Morus serrator",
    "Sula sula" = "Red-footed Booby|Sula sula",
    "Sula leucogaster" = "(B|b)rown (B|b)ooby|Sula leucogaster",
    "Sula nebouxii" = "(B|b)lue-footed (B|b)ooby|Sula nebouxii",
    "Sula variegata" = "Peruvian Booby|Sula variegata",
    "Sula dactylatra" = "Masked Booby|Sula dactylatra",
    "Sula granti" = "Nazca Booby|Sula granti",
    "Microcarbo coronatus" = "Crowned[- _](c|C)ormorant|Microcarbo coronatus",
    "Microcarbo africanus" = "Long-tailed[- _](c|C)ormorant|Microcarbo africanus",
    "Microcarbo pygmaeus" = "Pygmy[- _](c|C)ormorant|Microcarbo pygmaeus",
    "Microcarbo niger" = "Little[- _](c|C)ormorant|Microcarbo niger",
    "Microcarbo melanoleucos" = "Little Pied[- _](c|C)ormorant|Microcarbo melanoleucos",
    "Poikilocarbo gaimardi" = "Red-legged[- _](c|C)ormorant|Poikilocarbo gaimardi",
    "Leucocarbo magellanicus" = "Rock Shag|Leucocarbo magellanicus",
    "Leucocarbo bougainvilliorum" = "(G|g)uanay[ -_](c|C)ormorant|Leucocarbo bougainvilliorum",
    "Leucocarbo atriceps" = "Imperial Shag|Leucocarbo atriceps",
    "Leucocarbo verrucosus" = "Kerguelen Islands Shag|Leucocarbo verrucosus",
    "Leucocarbo carunculatus" = "Rough-faced Shag|Leucocarbo carunculatus",
    "Leucocarbo chalconotus" = "Stewart Island Shag|Leucocarbo chalconotus",
    "Leucocarbo onslowi" = "Chatham Islands Shag|Leucocarbo onslowi",
    "Leucocarbo campbelli" = "Campbell Island Shag|Leucocarbo campbelli",
    "Leucocarbo ranfurlyi" = "(B|b)ounty Islands Shag|Leucocarbo ranfurlyi",
    "Leucocarbo colensoi" = "Auckland Islands Shag|Leucocarbo colensoi",
    "Nannopterum auritum" = "(D|d)ouble[- _](c|C)rested[- _](c|C)ormorant|Nannopterum auritum",
    "Nannopterum brasilianum" = "(N|n)eotropic(al){0,1}[- _](c|C)ormorant|Nannopterum brasilianum",
    "Nannopterum harrisi" = "Flightless[- _](c|C)ormorant|Nannopterum harrisi",
    "Urile penicillatus" = "(B|b)randt('{0,1}s){0,1}[- _](c|C)ormorant|Urile penicillatus",
    "Urile pelagicus" = "Pelagic[- _](c|C)ormorant|Urile pelagicus",
    "Urile urile" = "Red-faced[- _](c|C)ormorant|Urile urile",
    "Urile perspicillatus" = "Spectacled[- _](c|C)ormorant|Urile perspicillatus",
    "Gulosus aristotelis" = "European Shag|Gulosus aristotelis",
    "Phalacrocorax carbo" = "((G|g)reat|White[- _]{0,1}breasted)[- _]{0,1}(C|c)ormorant|(P|p)halacrocorax[- _]{0,1}carbo",
    "Phalacrocorax capillatus" = "Japanese[- _](c|C)ormorant|(P|p)halacrocorax[- _]{0,1}capillatus",
    "Phalacrocorax capensis" = "(C|c)ape (C|c)ormorant|(P|p)halacrocorax[- _]{0,1}capensis",
    "Phalacrocorax nigrogularis" = "Socotra[- _](c|C)ormorant|(P|p)halacrocorax[- _]{0,1}nigrogularis",
    "Phalacrocorax neglectus" = "(B|b)ank[- _](c|C)ormorant|(P|p)halacrocorax[- _]{0,1}neglectus",
    "Phalacrocorax fuscicollis" = "Indian[- _](c|C)ormorant|(P|p)halacrocorax[- _]{0,1}fuscicollis",
    "Phalacrocorax sulcirostris" = "Little Black[- _](c|C)ormorant|(P|p)halacrocorax[- _]{0,1}sulcirostris",
    "Phalacrocorax fuscescens" = "(B|b)lack-faced[- _](c|C)ormorant|(P|p)halacrocorax[- _]{0,1}fuscescens",
    "Phalacrocorax varius" = "Great Pied[- _](c|C)ormorant|(P|p)halacrocorax[- _]{0,1}varius",
    "Phalacrocorax punctatus" = "Spotted Shag|(P|p)halacrocorax[- _]{0,1}punctatus",
    "Phalacrocorax featherstoni" = "Pitt Island Shag|(P|p)halacrocorax[- _]{0,1}featherstoni",
    'Sulidae sp.' = '^(g|G)annet{0,1}$',
    'Phalacrocoracidae sp.' = '^(c|C)ormorants{0,1}$'
    
  )
  
  z <- NA_character_
  
  for (i in 1:length(suliformes)) {
    if (any(grepl(suliformes[[i]], x))) {
      z <- names(suliformes)[[i]]
    }
  }
  
  return(z)
}


FormatPhaethontiformes <- function(x){
  phaethontiformes <- c(
    "Phaethon aethereus" = "Red-billed Tropicbird|Phaethon aethereus",
    "Phaethon rubricauda" = "Red-tailed Tropicbird|Phaethon rubricauda",
    "Phaethon lepturus" = "White-tailed Tropicbird|Phaethon lepturus"
  )
  z <- NA_character_
  for (i in 1:length(phaethontiformes)) {
    if (any(grepl(phaethontiformes[[i]], x))) {
      z <- names(phaethontiformes)[[i]]
    }
  }
  
  return(z)
}


#################################### END #######################################
################################################################################