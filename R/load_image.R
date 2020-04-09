load_image <- function(file, folder = "https://raw.githubusercontent.com/NUS-VIP/predicting-human-gaze-beyond-pixels/master/data/stimuli"){
  require(imager)
  
  img <- imager::load.image(file.path(folder, file))
  
  return(img)
}