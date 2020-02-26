html.vec <- readLines("https://lftdi.camden.rutgers.edu/provedit/files/")
zip.vec <- grep("zip", html.vec, value=TRUE)
zip.dt <- nc::capture_first_vec(
  zip.vec,
  nc::field("href", '="', '[^"]+'),
  '">',
  name=".*?",
  "<")
dir.create("zip")

for(zip.i in 1:nrow(zip.dt)){
  z <- zip.dt[zip.i]
  path <- file.path("zip", z$name)
  if(!file.exists(path)){
    cat(sprintf(
      "%4d / %4d %s -> %s\n",
      zip.i, nrow(zip.dt),
      z$href, path))
    download.file(URLencode(z$href), path)
  }
  unzip(path)
}

system("mkdir -p data && mv PROVEDIt_* data")
