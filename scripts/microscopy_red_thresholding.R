# install.packages('imager')
library(imager)

img.path <- "input/microscopy/221220/JPEGs/"
img.out <- "output/microscopy/221220/"
img.files <- list.files(path = img.path, pattern = '.jpg', all.files = T)


for (i in img.files) {
  img <- load.image(sprintf('%s%s',img.path, i))
  cn <- imsplit(img,"c")
  save.image(im = cn$`c = 1`, file = sprintf('%s%s',img.out, i))
}


