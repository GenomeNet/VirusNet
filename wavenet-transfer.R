library(deepG)
library(keras)
library(tensorflow)

tf$compat$v1$disable_eager_execution()

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

with(tf$device("/cpu:0"), {
  model <- create_model_wavenet(filters = 30, kernel_size = 16, residual_blocks = 2 ^ rep(0:9, 4), maxlen = 10000,
                                  input_tensor = NULL, initial_kernel_size = 3, initial_filters = 3 ^ 4,
                                  output_channels = 4, output_activation = "softmax", solver = "adam",
                                  learning.rate = 0.001, compile = F)
})

model <- multi_gpu_model(model, gpus = 2, cpu_merge = F, cpu_relocation = F)
model <- model %>% keras::compile(loss = "categorical_crossentropy",
                           optimizer = "adam", metrics = c("acc"))

keras::load_model_weights_hdf5(model, "~/workspace/checkpoints/Ep.006-val_loss1.29-val_acc0.359.hdf5")


deepG::trainNetwork(train_type = "lm",
             model = model,
             path = "~/workspace/training_data/bacteria_0_1/train",
             path.val = "~/workspace/training_data/bacteria_0_1/validation/",
             checkpoint_path = "~/workspace/checkpoints",
             validation.split = 0.1,
             run.name = "wave10k_new",
             batch.size = 40,
             epochs = 2,
             percentage_per_file = 0.33,
             max.queue.size = 500,
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
                           log = FALSE,
                           serialize_model = FALSE,
                           full_model = FALSE
             ),
             reverseComplements = FALSE,
       wavenet_format = T)


pretrained <- deepG::create_model_wavenet(filters = 30, kernel_size = 16, residual_blocks = 2 ^ rep(0:9, 4), maxlen = 10000,
                                  input_tensor = NULL, initial_kernel_size = 3, initial_filters = 3 ^ 4,
                                  output_channels = 4, output_activation = "softmax", solver = "adam",
                                  learning.rate = 0.001)

keras::load_model_weights_hdf5(model2, "~/workspace/checkpoints/Ep.006-val_loss1.29-val_acc0.359.hdf5")

keras::freeze_weights(model3)

base <- get_layer(pretrained, index = -1)$output
target <- base %>%
#layer_reshape(target_shape = c(1000, 40)) %>%
layer_dense(units = 100, trainable = T) %>%
   layer_activation("relu") %>%
   layer_dense(2, trainable = T) %>%
   layer_activation("softmax", trainable = T)

transfer_model <- keras_model(inputs = pretrained$input, output = target)

summary(transfer_model)

transfer_model %>% keras::compile(loss = "categorical_crossentropy",
                           optimizer = "adam", metrics = c("acc"))


deepG::trainNetwork(train_type = "label_folder",
             model = transfer_model,
             path = c("~/workspace/training_data/bacteria_0_1/train",
             "~/workspace/training_data/viral_0_1/train"),
             path.val = c("~/workspace/training_data/bacteria_0_1/validation/",
             "~/workspace/training_data/viral_0_1/validation"),
             checkpoint_path = "/workspace/checkpoints",
             validation.split = 0.2,
             run.name = "transfer_wave",
             batch.size = 10,
            percentage_per_file = 0.1,
             epochs = 10,
             max.queue.size = 50,
             lr.plateau.factor = 0.9,
             patience = 5,
             cooldown = 2,
             steps.per.epoch = 1000,
             step = 400,
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
             labelVocabulary = c("bacteria", "virus"),
             reverseComplements = FALSE,
       wavenet_format = F)
