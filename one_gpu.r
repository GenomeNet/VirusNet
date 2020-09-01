library(keras)
library(tensorflow)
library(deepG)
library(wavenet)

tf$config$list_physical_devices("GPU")
#' create wavenet model
#'  
#' @inheritParams wavenet::wavenet
#' @export
create_model_wavenet <- function(filters = 16, kernel_size = 2, residual_blocks, maxlen,
                                 input_tensor = NULL, initial_kernel_size = 32, initial_filters = 32,
                                 output_channels = 4, output_activation = "softmax", solver = "adam",
                                 learning.rate = 1e-3, compile = T) {

  model <- wavenet::wavenet(filters = filters, kernel_size = kernel_size, residual_blocks = residual_blocks,
                            input_shape = list(maxlen, output_channels), input_tensor = input_tensor, initial_kernel_size = initial_kernel_size,
                            initial_filters = initial_filters, output_channels = output_channels, output_activation = "softmax")
  if (solver == "adam") {
    optimizer <- keras::optimizer_adam(lr = learning.rate)
  }
  if (solver == "adagrad") {
    optimizer <- keras::optimizer_adagrad(lr = learning.rate)
  }
  if (solver == "rmsprop") {
    optimizer <- keras::optimizer_rmsprop(lr = learning.rate)
  }
  if (solver == "sgd") {
    optimizer <- keras::optimizer_sgd(lr = learning.rate)
  }
  if (compile)
    model %>% keras::compile(loss = "categorical_crossentropy",
                           optimizer = optimizer, metrics = c("acc"))

  argg <- c(as.list(environment()))
  argg["model"] <- NULL
  argg["optimizer"] <- NULL
  argg["residual_blocks"] <- paste(as.character(residual_blocks), collapse = " ")
  model$hparam <- argg
  model
}


model <- create_model_wavenet(filters = 30, kernel_size = 16, residual_blocks = 2 ^ rep(0:9, 4), maxlen = 1000,
                                  input_tensor = NULL, initial_kernel_size = 3, initial_filters = 3 ^ 4,
                                  output_channels = 4, output_activation = "softmax", solver = "adam",
                                  learning.rate = 0.001, compile = F)

model <- model %>% keras::compile(loss = "categorical_crossentropy",
                           optimizer = "adam", metrics = c("acc"))

#keras::load_model_weights_hdf5(model, "~/workspace/checkpoints/Ep.006-val_loss1.29-val_acc0.359.hdf5")

deepG::trainNetwork(train_type = "lm",
             model = model,
             path = "~/workspace/training_data/bacteria_0_1/train",
             path.val = "~/workspace/training_data/bacteria_0_1/validation/",
             checkpoint_path = "~/workspace/checkpoints",
             validation.split = 0.1,
             run.name = "test2",
             batch.size = 500,
             epochs = 2,
             percentage_per_file = 0.33,
             max.queue.size = 1000,
             lr.plateau.factor = 0.9,
             patience = 5,
             cooldown = 2,
             steps.per.epoch = 1000,
             step = 5000,
             vocabulary = c("A", "C", "T", "G"),
             tensorboard.log = "/workspace/tensorboard",
             shuffleFastaEntries = TRUE,
             output = list(none = FALSE,
                           checkpoints = TRUE,
                           tensorboard = FALSE,
                           log = TRUE,
                           serialize_model = FALSE,
                           full_model = FALSE
             ),
             reverseComplements = TRUE,
       wavenet_format = T)
