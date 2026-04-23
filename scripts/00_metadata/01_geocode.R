# scripts/00_metadata/01_geocode.R
#
# PURPOSE: Geocode unique GISAID location strings using Nominatim (OpenStreetMap)
#          via the tidygeocoder package. Results are cached so re-runs only query
#          new locations. Coordinates are written back into all gene metadata CSVs.
#
# URBAN FALLBACK:
#   Location strings with only continent/country (2 components) or
#   continent/country/state (3 components) get geocoded to geographic centroids,
#   which are typically rural and have near-zero population density. These are
#   replaced with coordinates for the largest city in that state/country via a
#   lookup table. Entries not in the lookup table fall back to a Nominatim
#   query with city-type preference.
#
# INPUT:   data/metadata/{gene}_metadata.csv  (run 00_parse_gisaid.R first)
# OUTPUT:  data/raw/location_coords_cache.csv — persistent geocoding cache
#          data/raw/geocode_failures.txt       — locations that could not be resolved
#          data/metadata/{gene}_metadata.csv   — lat/lon columns filled in
#
# RATE LIMIT: Nominatim requires <= 1 request/sec. tidygeocoder enforces this
#             automatically. For ~1000 unique locations expect ~20 min.
#
# USAGE:
#   Rscript scripts/00_metadata/01_geocode.R

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidygeocoder)
})

source(here("config.R"))

ALL_GENES    <- c(GENES, "PA_H1N1", "PA_H3N2")
cache_path   <- here("data", "raw", "location_coords_cache.csv")
failures_path <- here("data", "raw", "geocode_failures.txt")

max_fallback_levels <- 3   # how many location components to strip before giving up

# ── Urban fallback lookup tables ──────────────────────────────────────────────
# For location strings with only continent/country/state (3 components) or
# continent/country (2 components), Nominatim returns a geographic centroid
# that is typically rural. These tables map (country, admin1) → the Nominatim
# query string for the largest/most representative city in that region.
#
# country names match GISAID location string format (e.g. "United States",
# not "USA"). admin1 matching is case-insensitive and strips common suffixes
# ("Province", "SAR", "Municipality", "Autonomous Region", etc.).

URBAN_STATE_LOOKUP <- tribble(
  ~country,           ~admin1,                              ~city_query,
  # United States — largest city by state
  "United States", "Alabama",                "Birmingham, Alabama, United States",
  "United States", "Alaska",                 "Anchorage, Alaska, United States",
  "United States", "Arizona",                "Phoenix, Arizona, United States",
  "United States", "Arkansas",               "Little Rock, Arkansas, United States",
  "United States", "California",             "Los Angeles, California, United States",
  "United States", "Colorado",               "Denver, Colorado, United States",
  "United States", "Connecticut",            "Bridgeport, Connecticut, United States",
  "United States", "Delaware",               "Wilmington, Delaware, United States",
  "United States", "District of Columbia",   "Washington DC, United States",
  "United States", "Florida",                "Jacksonville, Florida, United States",
  "United States", "Georgia",                "Atlanta, Georgia, United States",
  "United States", "Hawaii",                 "Honolulu, Hawaii, United States",
  "United States", "Idaho",                  "Boise, Idaho, United States",
  "United States", "Illinois",               "Chicago, Illinois, United States",
  "United States", "Indiana",                "Indianapolis, Indiana, United States",
  "United States", "Iowa",                   "Des Moines, Iowa, United States",
  "United States", "Kansas",                 "Wichita, Kansas, United States",
  "United States", "Kentucky",               "Louisville, Kentucky, United States",
  "United States", "Louisiana",              "New Orleans, Louisiana, United States",
  "United States", "Maine",                  "Portland, Maine, United States",
  "United States", "Maryland",               "Baltimore, Maryland, United States",
  "United States", "Massachusetts",          "Boston, Massachusetts, United States",
  "United States", "Michigan",               "Detroit, Michigan, United States",
  "United States", "Minnesota",              "Minneapolis, Minnesota, United States",
  "United States", "Mississippi",            "Jackson, Mississippi, United States",
  "United States", "Missouri",               "Kansas City, Missouri, United States",
  "United States", "Montana",                "Billings, Montana, United States",
  "United States", "Nebraska",               "Omaha, Nebraska, United States",
  "United States", "Nevada",                 "Las Vegas, Nevada, United States",
  "United States", "New Hampshire",          "Manchester, New Hampshire, United States",
  "United States", "New Jersey",             "Newark, New Jersey, United States",
  "United States", "New Mexico",             "Albuquerque, New Mexico, United States",
  "United States", "New York",               "New York City, United States",
  "United States", "North Carolina",         "Charlotte, North Carolina, United States",
  "United States", "North Dakota",           "Fargo, North Dakota, United States",
  "United States", "Ohio",                   "Columbus, Ohio, United States",
  "United States", "Oklahoma",               "Oklahoma City, Oklahoma, United States",
  "United States", "Oregon",                 "Portland, Oregon, United States",
  "United States", "Pennsylvania",           "Philadelphia, Pennsylvania, United States",
  "United States", "Rhode Island",           "Providence, Rhode Island, United States",
  "United States", "South Carolina",         "Columbia, South Carolina, United States",
  "United States", "South Dakota",           "Sioux Falls, South Dakota, United States",
  "United States", "Tennessee",              "Memphis, Tennessee, United States",
  "United States", "Texas",                  "Houston, Texas, United States",
  "United States", "Utah",                   "Salt Lake City, Utah, United States",
  "United States", "Vermont",                "Burlington, Vermont, United States",
  "United States", "Virginia",               "Virginia Beach, Virginia, United States",
  "United States", "Washington",             "Seattle, Washington, United States",
  "United States", "West Virginia",          "Charleston, West Virginia, United States",
  "United States", "Wisconsin",              "Milwaukee, Wisconsin, United States",
  "United States", "Wyoming",                "Cheyenne, Wyoming, United States",
  # China — provincial capital or largest city; admin1 matched after suffix-stripping
  "China", "Beijing",            "Beijing, China",
  "China", "Tianjin",            "Tianjin, China",
  "China", "Shanghai",           "Shanghai, China",
  "China", "Chongqing",          "Chongqing, China",
  "China", "Hong Kong",          "Hong Kong",
  "China", "Hebei",              "Shijiazhuang, Hebei, China",
  "China", "Shanxi",             "Taiyuan, Shanxi, China",
  "China", "Inner Mongolia",     "Hohhot, Inner Mongolia, China",
  "China", "Liaoning",           "Shenyang, Liaoning, China",
  "China", "Jilin",              "Changchun, Jilin, China",
  "China", "Heilongjiang",       "Harbin, Heilongjiang, China",
  "China", "Jiangsu",            "Nanjing, Jiangsu, China",
  "China", "Zhejiang",           "Hangzhou, Zhejiang, China",
  "China", "Anhui",              "Hefei, Anhui, China",
  "China", "Fujian",             "Fuzhou, Fujian, China",
  "China", "Jiangxi",            "Nanchang, Jiangxi, China",
  "China", "Shandong",           "Jinan, Shandong, China",
  "China", "Henan",              "Zhengzhou, Henan, China",
  "China", "Hubei",              "Wuhan, Hubei, China",
  "China", "Hunan",              "Changsha, Hunan, China",
  "China", "Guangdong",          "Guangzhou, Guangdong, China",
  "China", "Guangxi",            "Nanning, Guangxi, China",
  "China", "Hainan",             "Haikou, Hainan, China",
  "China", "Sichuan",            "Chengdu, Sichuan, China",
  "China", "Guizhou",            "Guiyang, Guizhou, China",
  "China", "Yunnan",             "Kunming, Yunnan, China",
  "China", "Tibet",              "Lhasa, Tibet, China",
  "China", "Shaanxi",            "Xi'an, Shaanxi, China",
  "China", "Gansu",              "Lanzhou, Gansu, China",
  "China", "Qinghai",            "Xining, Qinghai, China",
  "China", "Ningxia",            "Yinchuan, Ningxia, China",
  "China", "Xinjiang",           "Urumqi, Xinjiang, China",
  # Australia — largest city by state/territory
  "Australia", "New South Wales",             "Sydney, Australia",
  "Australia", "Victoria",                    "Melbourne, Australia",
  "Australia", "Queensland",                  "Brisbane, Australia",
  "Australia", "Western Australia",           "Perth, Australia",
  "Australia", "South Australia",             "Adelaide, Australia",
  "Australia", "Tasmania",                    "Hobart, Australia",
  "Australia", "Australian Capital Territory","Canberra, Australia",
  "Australia", "Northern Territory",          "Darwin, Australia",
  # Canada — largest city by province
  "Canada", "Ontario",                    "Toronto, Canada",
  "Canada", "Quebec",                     "Montreal, Canada",
  "Canada", "British Columbia",           "Vancouver, Canada",
  "Canada", "Alberta",                    "Calgary, Canada",
  "Canada", "Manitoba",                   "Winnipeg, Canada",
  "Canada", "Saskatchewan",              "Saskatoon, Canada",
  "Canada", "Nova Scotia",               "Halifax, Canada",
  "Canada", "New Brunswick",             "Moncton, Canada",
  "Canada", "Newfoundland and Labrador", "St. John's, Newfoundland, Canada",
  "Canada", "Prince Edward Island",      "Charlottetown, Canada",
  "Canada", "Northwest Territories",     "Yellowknife, Canada",
  "Canada", "Yukon",                     "Whitehorse, Canada",
  "Canada", "Nunavut",                   "Iqaluit, Canada",
  # Brazil — largest city by state
  "Brazil", "Sao Paulo",           "Sao Paulo, Brazil",
  "Brazil", "Rio de Janeiro",      "Rio de Janeiro, Brazil",
  "Brazil", "Minas Gerais",        "Belo Horizonte, Brazil",
  "Brazil", "Bahia",               "Salvador, Brazil",
  "Brazil", "Rio Grande do Sul",   "Porto Alegre, Brazil",
  "Brazil", "Parana",              "Curitiba, Brazil",
  "Brazil", "Pernambuco",          "Recife, Brazil",
  "Brazil", "Ceara",               "Fortaleza, Brazil",
  "Brazil", "Goias",               "Goiania, Brazil",
  "Brazil", "Amazonas",            "Manaus, Brazil",
  "Brazil", "Para",                "Belem, Brazil",
  "Brazil", "Maranhao",            "Sao Luis, Brazil",
  "Brazil", "Santa Catarina",      "Joinville, Brazil",
  "Brazil", "Espirito Santo",      "Vitoria, Brazil",
  "Brazil", "Mato Grosso",         "Cuiaba, Brazil",
  "Brazil", "Mato Grosso do Sul",  "Campo Grande, Brazil",
  # India — largest city by state
  "India", "Maharashtra",                  "Mumbai, India",
  "India", "Karnataka",                    "Bengaluru, India",
  "India", "Tamil Nadu",                   "Chennai, India",
  "India", "West Bengal",                  "Kolkata, India",
  "India", "Delhi",                        "New Delhi, India",
  "India", "National Capital Territory of Delhi", "New Delhi, India",
  "India", "Telangana",                    "Hyderabad, India",
  "India", "Gujarat",                      "Ahmedabad, India",
  "India", "Uttar Pradesh",                "Lucknow, India",
  "India", "Rajasthan",                    "Jaipur, India",
  "India", "Punjab",                       "Ludhiana, India",
  "India", "Andhra Pradesh",               "Visakhapatnam, India",
  "India", "Madhya Pradesh",               "Indore, India",
  "India", "Bihar",                        "Patna, India",
  "India", "Odisha",                       "Bhubaneswar, India",
  "India", "Kerala",                       "Thiruvananthapuram, India",
  "India", "Haryana",                      "Faridabad, India",
  "India", "Uttarakhand",                  "Dehradun, India",
  # United Kingdom — by nation/region
  "United Kingdom", "England",          "London, United Kingdom",
  "United Kingdom", "Scotland",         "Glasgow, United Kingdom",
  "United Kingdom", "Wales",            "Cardiff, United Kingdom",
  "United Kingdom", "Northern Ireland", "Belfast, United Kingdom",
  # Germany — by state (Bundesland)
  "Germany", "Baden-Wurttemberg",  "Stuttgart, Germany",
  "Germany", "Bavaria",            "Munich, Germany",
  "Germany", "Berlin",             "Berlin, Germany",
  "Germany", "Brandenburg",        "Potsdam, Germany",
  "Germany", "Bremen",             "Bremen, Germany",
  "Germany", "Hamburg",            "Hamburg, Germany",
  "Germany", "Hesse",              "Frankfurt, Germany",
  "Germany", "Lower Saxony",       "Hanover, Germany",
  "Germany", "Mecklenburg-Vorpommern", "Rostock, Germany",
  "Germany", "North Rhine-Westphalia", "Cologne, Germany",
  "Germany", "Rhineland-Palatinate",  "Mainz, Germany",
  "Germany", "Saarland",           "Saarbrucken, Germany",
  "Germany", "Saxony",             "Leipzig, Germany",
  "Germany", "Saxony-Anhalt",      "Halle, Germany",
  "Germany", "Schleswig-Holstein", "Kiel, Germany",
  "Germany", "Thuringia",          "Erfurt, Germany",
  # South Africa — by province
  "South Africa", "Gauteng",                  "Johannesburg, South Africa",
  "South Africa", "KwaZulu Natal",            "Durban, South Africa",
  "South Africa", "Western Cape",             "Cape Town, South Africa",
  "South Africa", "Province of the Western Cape", "Cape Town, South Africa",
  "South Africa", "Eastern Cape",             "East London, South Africa",
  "South Africa", "Mpumalanga",               "Mbombela, South Africa",
  "South Africa", "Limpopo",                  "Polokwane, South Africa",
  "South Africa", "North West",               "Rustenburg, South Africa",
  "South Africa", "Free State",               "Bloemfontein, South Africa",
  "South Africa", "Northern Cape",            "Kimberley, South Africa",
  # Thailand — largest city per province (province name = city for major ones)
  "Thailand", "Bangkok",       "Bangkok, Thailand",
  "Thailand", "Chiang Mai",    "Chiang Mai, Thailand",
  "Thailand", "Chiang Rai",    "Chiang Rai, Thailand",
  "Thailand", "Surat Thani",   "Surat Thani, Thailand",
  "Thailand", "Nakhon Phanom", "Nakhon Phanom, Thailand",
  "Thailand", "Chanthaburi",   "Chanthaburi, Thailand",
  "Thailand", "Yasothon",      "Yasothon, Thailand",
  "Thailand", "Satun",         "Satun, Thailand",
  # Japan — prefecture name is the city (most cases)
  "Japan", "Chiba",     "Chiba, Japan",
  "Japan", "Gunma",     "Maebashi, Japan",
  "Japan", "Iwate",     "Morioka, Japan",
  "Japan", "Kyoto",     "Kyoto, Japan",
  "Japan", "Nagasaki",  "Nagasaki, Japan",
  "Japan", "Niigata",   "Niigata, Japan",
  "Japan", "Saitama",   "Saitama, Japan",
  # Russian Federation — by region
  "Russian Federation", "Stavropol Krai",   "Stavropol, Russia",
  "Russian Federation", "Respublica Krym",  "Simferopol, Russia",
  # Kenya — by province (assign to nearest major city)
  "Kenya", "Central Province",  "Nyeri, Kenya",
  "Kenya", "Coast Province",    "Mombasa, Kenya",
  "Kenya", "Nyanza Province",   "Kisumu, Kenya",
  "Kenya", "Western Province",  "Kakamega, Kenya",
  "Kenya", "Nairobi",           "Nairobi, Kenya",
  "Kenya", "Nakuru",            "Nakuru, Kenya",
  "Kenya", "Nyeri",             "Nyeri, Kenya",
  "Kenya", "Kilifi",            "Kilifi, Kenya"
)

# Country-level fallback (2-component locations: continent/country)
# Maps country name → largest/capital city query
COUNTRY_CAPITALS <- tribble(
  ~country,                    ~city_query,
  "United States",             "New York City, United States",
  "China",                     "Beijing, China",
  "Russia",                    "Moscow, Russia",
  "Russian Federation",        "Moscow, Russia",
  "United Kingdom",            "London, United Kingdom",
  "Germany",                   "Berlin, Germany",
  "France",                    "Paris, France",
  "Italy",                     "Rome, Italy",
  "Spain",                     "Madrid, Spain",
  "Brazil",                    "Sao Paulo, Brazil",
  "India",                     "Mumbai, India",
  "Australia",                 "Sydney, Australia",
  "Canada",                    "Toronto, Canada",
  "Japan",                     "Tokyo, Japan",
  "South Korea",               "Seoul, South Korea",
  "Mexico",                    "Mexico City, Mexico",
  "Indonesia",                 "Jakarta, Indonesia",
  "Turkey",                    "Istanbul, Turkey",
  "Netherlands",               "Amsterdam, Netherlands",
  "Belgium",                   "Brussels, Belgium",
  "Sweden",                    "Stockholm, Sweden",
  "Norway",                    "Oslo, Norway",
  "Denmark",                   "Copenhagen, Denmark",
  "Finland",                   "Helsinki, Finland",
  "Poland",                    "Warsaw, Poland",
  "Czech Republic",            "Prague, Czech Republic",
  "Austria",                   "Vienna, Austria",
  "Switzerland",               "Zurich, Switzerland",
  "Portugal",                  "Lisbon, Portugal",
  "Greece",                    "Athens, Greece",
  "Romania",                   "Bucharest, Romania",
  "Hungary",                   "Budapest, Hungary",
  "Ukraine",                   "Kyiv, Ukraine",
  "Argentina",                 "Buenos Aires, Argentina",
  "Colombia",                  "Bogota, Colombia",
  "Chile",                     "Santiago, Chile",
  "Peru",                      "Lima, Peru",
  "Venezuela",                 "Caracas, Venezuela",
  "South Africa",              "Johannesburg, South Africa",
  "Nigeria",                   "Lagos, Nigeria",
  "Egypt",                     "Cairo, Egypt",
  "Kenya",                     "Nairobi, Kenya",
  "Ethiopia",                  "Addis Ababa, Ethiopia",
  "Ghana",                     "Accra, Ghana",
  "Tanzania",                  "Dar es Salaam, Tanzania",
  "Morocco",                   "Casablanca, Morocco",
  "Tunisia",                   "Tunis, Tunisia",
  "Thailand",                  "Bangkok, Thailand",
  "Vietnam",                   "Ho Chi Minh City, Vietnam",
  "Malaysia",                  "Kuala Lumpur, Malaysia",
  "Philippines",               "Manila, Philippines",
  "Singapore",                 "Singapore",
  "Bangladesh",                "Dhaka, Bangladesh",
  "Pakistan",                  "Karachi, Pakistan",
  "Iran",                      "Tehran, Iran",
  "Saudi Arabia",              "Riyadh, Saudi Arabia",
  "Israel",                    "Tel Aviv, Israel",
  "New Zealand",               "Auckland, New Zealand",
  "Taiwan",                    "Taipei, Taiwan",
  "Hong Kong",                 "Hong Kong",
  "Serbia",                    "Belgrade, Serbia",
  "Croatia",                   "Zagreb, Croatia",
  "Slovakia",                  "Bratislava, Slovakia",
  "Bulgaria",                  "Sofia, Bulgaria",
  "Ireland",                   "Dublin, Ireland",
  "Iceland",                   "Reykjavik, Iceland",
  "Luxembourg",                "Luxembourg City, Luxembourg",
  "Cuba",                      "Havana, Cuba",
  "Dominican Republic",        "Santo Domingo, Dominican Republic",
  "Guatemala",                 "Guatemala City, Guatemala",
  "Costa Rica",                "San Jose, Costa Rica",
  "Panama",                    "Panama City, Panama",
  "Ecuador",                   "Guayaquil, Ecuador",
  "Bolivia",                   "Santa Cruz, Bolivia",
  "Uruguay",                   "Montevideo, Uruguay",
  "Cameroon",                  "Douala, Cameroon",
  "Senegal",                   "Dakar, Senegal",
  "Zimbabwe",                  "Harare, Zimbabwe",
  "Zambia",                    "Lusaka, Zambia",
  "Uganda",                    "Kampala, Uganda",
  "Rwanda",                    "Kigali, Rwanda",
  "Myanmar",                   "Yangon, Myanmar",
  "Cambodia",                  "Phnom Penh, Cambodia",
  "Sri Lanka",                 "Colombo, Sri Lanka",
  "Nepal",                     "Kathmandu, Nepal",
  "Kazakhstan",                "Almaty, Kazakhstan",
  "Lebanon",                   "Beirut, Lebanon",
  # Africa — missing from original list
  "Georgia",                   "Tbilisi, Georgia",  # country, not US state
  "Tanzania, United Republic of", "Dar es Salaam, Tanzania",
  "Cote d'Ivoire",             "Abidjan, Cote d'Ivoire",
  "Burkina Faso",              "Ouagadougou, Burkina Faso",
  "Central African Republic",  "Bangui, Central African Republic",
  "Gambia",                    "Banjul, Gambia",
  "Madagascar",                "Antananarivo, Madagascar",
  "Mali",                      "Bamako, Mali",
  "Mauritius",                 "Port Louis, Mauritius",
  "Mozambique",                "Maputo, Mozambique",
  "Niger",                     "Niamey, Niger",
  "Sierra Leone",              "Freetown, Sierra Leone",
  "Sudan",                     "Khartoum, Sudan",
  "Togo",                      "Lome, Togo",
  # Asia — missing
  "Afghanistan",               "Kabul, Afghanistan",
  "Armenia",                   "Yerevan, Armenia",
  "Bhutan",                    "Thimphu, Bhutan",
  "Brunei",                    "Bandar Seri Begawan, Brunei",
  "Hong Kong (SAR)",           "Hong Kong",
  "Iraq",                      "Baghdad, Iraq",
  "Jordan",                    "Amman, Jordan",
  "Kuwait",                    "Kuwait City, Kuwait",
  "Kyrgyzstan",                "Bishkek, Kyrgyzstan",
  "Mongolia",                  "Ulaanbaatar, Mongolia",
  "Oman",                      "Muscat, Oman",
  "Palestinian Territory",     "Ramallah, Palestinian Territory",
  "Qatar",                     "Doha, Qatar",
  "Tajikistan",                "Dushanbe, Tajikistan",
  "Turkmenistan",              "Ashgabat, Turkmenistan",
  "United Arab Emirates",      "Dubai, United Arab Emirates",
  # Europe — missing
  "Albania",                   "Tirana, Albania",
  "Bosnia",                    "Sarajevo, Bosnia and Herzegovina",
  "Bosnia and Herzegovina",    "Sarajevo, Bosnia and Herzegovina",
  "Cyprus",                    "Nicosia, Cyprus",
  "Estonia",                   "Tallinn, Estonia",
  "Kosovo",                    "Pristina, Kosovo",
  "Latvia",                    "Riga, Latvia",
  "Lithuania",                 "Vilnius, Lithuania",
  "Montenegro",                "Podgorica, Montenegro",
  # Caribbean and Central America — missing
  "Aruba",                     "Oranjestad, Aruba",
  "Belize",                    "Belize City, Belize",
  "Dominica",                  "Roseau, Dominica",
  "Grenada",                   "Saint George's, Grenada",
  "Guadeloupe",                "Pointe-a-Pitre, Guadeloupe",
  "Haiti",                     "Port-au-Prince, Haiti",
  "Honduras",                  "Tegucigalpa, Honduras",
  "Jamaica",                   "Kingston, Jamaica",
  "Montserrat",                "Plymouth, Montserrat",
  "Nicaragua",                 "Managua, Nicaragua",
  "Puerto Rico",               "San Juan, Puerto Rico",
  "Saint Lucia",               "Castries, Saint Lucia",
  "Trinidad and Tobago",       "Port of Spain, Trinidad and Tobago",
  # Oceania — missing
  "New Caledonia",             "Noumea, New Caledonia",
  "Papua New Guinea",          "Port Moresby, Papua New Guinea",
  "Timor-Leste",               "Dili, Timor-Leste",
  # South America — missing
  "French Guiana",             "Cayenne, French Guiana",
  "Guyana",                    "Georgetown, Guyana",
  "Paraguay",                  "Asuncion, Paraguay",
  "Suriname",                  "Paramaribo, Suriname"
)

# ── Helper: parse GISAID location string ──────────────────────────────────────
# Returns a named list: continent, country, admin1, city (any may be NA).
parse_location <- function(location) {
  parts <- str_split(location, " / ")[[1]] |> str_trim()
  list(
    continent = if (length(parts) >= 1) parts[1] else NA_character_,
    country   = if (length(parts) >= 2) parts[2] else NA_character_,
    admin1    = if (length(parts) >= 3) parts[3] else NA_character_,
    city      = if (length(parts) >= 4) parts[4] else NA_character_,
    n_parts   = length(parts)
  )
}

# Strip common Chinese administrative suffixes for lookup matching.
normalise_admin1 <- function(x) {
  x |>
    str_trim() |>
    str_remove_all(regex(
      "\\s*(Province|SAR|Municipality|Autonomous Region|Uyghur Autonomous Region|Special Administrative Region)\\s*$",
      ignore_case = TRUE
    ))
}

# ── Helper: look up urban city query for a location string ────────────────────
# Returns the Nominatim query string to use, or NA if not in lookup tables.
lookup_urban_city <- function(location) {
  loc <- parse_location(location)

  if (is.na(loc$country)) return(NA_character_)

  if (!is.na(loc$admin1)) {
    # 3-component (or more with failed specific geocode): try state-level lookup
    admin1_norm <- normalise_admin1(loc$admin1)
    hit <- URBAN_STATE_LOOKUP |>
      filter(
        str_to_lower(country) == str_to_lower(loc$country),
        str_to_lower(admin1)  == str_to_lower(admin1_norm)
      )
    if (nrow(hit) > 0) return(hit$city_query[1])
  }

  # Country-level fallback
  hit <- COUNTRY_CAPITALS |>
    filter(str_to_lower(country) == str_to_lower(loc$country))
  if (nrow(hit) > 0) return(hit$city_query[1])

  NA_character_
}

# ── Step 1: collect unique location strings ───────────────────────────────────
message("Loading metadata and collecting unique locations...")

all_locations <- character(0)
for (g in ALL_GENES) {
  path <- file.path(PATHS$metadata, paste0(g, "_metadata.csv"))
  if (!file.exists(path)) {
    warning("Missing metadata, skipping: ", path)
    next
  }
  locs <- read_csv(path, col_select = "location",
                   show_col_types = FALSE)[["location"]]
  all_locations <- c(all_locations, locs)
}

unique_locs <- unique(all_locations[!is.na(all_locations) & all_locations != ""])
message("  Unique location strings: ", length(unique_locs))

# ── Step 2: load or initialise cache ─────────────────────────────────────────
if (file.exists(cache_path)) {
  cache <- read_csv(cache_path, show_col_types = FALSE)
  message("  Loaded cache: ", nrow(cache), " previously geocoded locations")
  # Add new columns if cache pre-dates them; fill any NAs left by old runs
  if (!"geocode_method" %in% names(cache))
    cache$geocode_method <- "nominatim"
  cache$geocode_method[is.na(cache$geocode_method)] <- "nominatim"
  if (!"urban_city" %in% names(cache))
    cache$urban_city <- NA_character_
} else {
  cache <- tibble(
    location      = character(),
    lat           = numeric(),
    lon           = numeric(),
    geocode_level = integer(),
    query_used    = character(),
    geocode_method = character(),
    urban_city    = character()
  )
  message("  No cache found — starting fresh")
}

to_geocode <- setdiff(unique_locs, cache$location)
message("  Locations to geocode: ", length(to_geocode))

# ── Step 3: geocode with fallback ─────────────────────────────────────────────
build_query <- function(location, drop_n = 0) {
  parts <- str_split(location, " / ")[[1]] |> str_trim()
  parts <- parts[-1]                          # drop continent
  if (drop_n > 0) parts <- head(parts, -drop_n)
  if (length(parts) == 0) return(NA_character_)
  paste(rev(parts), collapse = ", ")
}

geocode_one <- function(location) {
  for (level in 0:max_fallback_levels) {
    query <- build_query(location, drop_n = level)
    if (is.na(query)) break

    result <- geo(address = query, method = "osm",
                  lat = "lat", long = "lon",
                  quiet = TRUE)

    if (!is.na(result$lat)) {
      return(tibble(
        location      = location,
        lat           = result$lat,
        lon           = result$lon,
        geocode_level = as.integer(level),
        query_used    = query,
        geocode_method = "nominatim",
        urban_city    = NA_character_
      ))
    }
  }

  return(tibble(
    location      = location,
    lat           = NA_real_,
    lon           = NA_real_,
    geocode_level = NA_integer_,
    query_used    = NA_character_,
    geocode_method = "nominatim",
    urban_city    = NA_character_
  ))
}

if (length(to_geocode) > 0) {
  message("\nGeocoding ", length(to_geocode), " locations via Nominatim ",
          "(~1 req/sec, this may take a while)...")

  new_results <- vector("list", length(to_geocode))
  for (i in seq_along(to_geocode)) {
    new_results[[i]] <- geocode_one(to_geocode[i])
    if (i %% 50 == 0) {
      message("  ", i, " / ", length(to_geocode), " done")
    }
  }
  new_results <- bind_rows(new_results)

  n_ok   <- sum(!is.na(new_results$lat))
  n_fail <- sum(is.na(new_results$lat))
  message("  Resolved: ", n_ok, "  Failed: ", n_fail)

  cache <- bind_rows(cache, new_results)
  write_csv(cache, cache_path)
  message("  Cache updated: ", nrow(cache), " total entries → ", cache_path)

  failures <- new_results$location[is.na(new_results$lat)]
  if (length(failures) > 0) {
    writeLines(failures, failures_path)
    message("  Failures logged to: ", failures_path)
  }
} else {
  message("  All locations already cached — nothing to query.")
}

# ── Step 4b: urban fallback for imprecise locations ───────────────────────────
# Locations with only continent/country (2 parts) or continent/country/state
# (3 parts) resolve to geographic centroids, which are typically rural.
# Replace these with coordinates for the largest city in that state/country.
# Also applied to 4+-component locations that fell back to state level or worse.
message("\nApplying urban fallback for imprecise locations...")

needs_fallback <- cache |>
  filter(!is.na(lat)) |>
  rowwise() |>
  mutate(n_parts = length(str_split(location, " / ")[[1]])) |>
  ungroup() |>
  filter(
    geocode_method != "urban_fallback",
    (n_parts <= 3) | (n_parts >= 4 & geocode_level >= 2)
  ) |>
  filter(!is.na(location))

message("  ", nrow(needs_fallback), " locations eligible for urban fallback")

if (nrow(needs_fallback) > 0) {
  needs_fallback <- needs_fallback |>
    rowwise() |>
    mutate(city_query = lookup_urban_city(location)) |>
    ungroup() |>
    filter(!is.na(city_query))

  message("  ", nrow(needs_fallback), " locations found in urban lookup table")

  # Geocode each unique city query (deduplicated to save API calls)
  unique_queries <- unique(needs_fallback$city_query)
  city_coords    <- vector("list", length(unique_queries))

  message("  Geocoding ", length(unique_queries),
          " unique urban centre queries via Nominatim...")

  for (i in seq_along(unique_queries)) {
    q      <- unique_queries[i]
    result <- geo(address = q, method = "osm",
                  lat = "lat", long = "lon", quiet = TRUE)
    city_coords[[i]] <- tibble(city_query = q,
                               city_lat   = result$lat,
                               city_lon   = result$lon)
    if (i %% 20 == 0)
      message("  ", i, " / ", length(unique_queries), " city queries done")
  }
  city_coords <- bind_rows(city_coords)

  n_city_ok <- sum(!is.na(city_coords$city_lat))
  message("  Urban city geocodes resolved: ", n_city_ok, " / ",
          nrow(city_coords))

  # Join city coordinates back and update the cache
  fallback_resolved <- needs_fallback |>
    left_join(city_coords, by = "city_query") |>
    filter(!is.na(city_lat))

  for (i in seq_len(nrow(fallback_resolved))) {
    loc <- fallback_resolved$location[i]
    cache <- cache |>
      rows_update(
        tibble(
          location       = loc,
          lat            = fallback_resolved$city_lat[i],
          lon            = fallback_resolved$city_lon[i],
          geocode_method = "urban_fallback",
          urban_city     = fallback_resolved$city_query[i]
        ),
        by = "location",
        unmatched = "ignore"
      )
  }

  write_csv(cache, cache_path)
  message("  Cache updated with ", nrow(fallback_resolved),
          " urban fallback coordinates → ", cache_path)

  n_no_lookup <- nrow(needs_fallback |>
    anti_join(fallback_resolved, by = "location"))
  if (n_no_lookup > 0)
    message("  NOTE: ", n_no_lookup,
            " locations had no lookup entry — kept original centroid coords.")
}

# ── Step 5: join coords back to all gene metadata CSVs ────────────────────────
message("\nWriting coordinates back to metadata CSVs...")

coords_lookup <- cache |>
  select(location, lat, lon) |>
  filter(!is.na(lat))

for (g in ALL_GENES) {
  path <- file.path(PATHS$metadata, paste0(g, "_metadata.csv"))
  if (!file.exists(path)) next

  df <- read_csv(path, show_col_types = FALSE)

  df <- df |>
    select(-latitude, -longitude) |>
    left_join(coords_lookup, by = "location") |>
    rename(latitude = lat, longitude = lon)

  df <- df |>
    select(accession, strain, subtype, lineage, passage_history,
           collection_date, year, month,
           location, country, region,
           latitude, longitude,
           host)

  write_csv(df, path)
  n_coords <- sum(!is.na(df$latitude))
  message("  [", g, "] ", n_coords, " / ", nrow(df),
          " rows with coordinates → ", basename(path))
}

message("\nDone.")
