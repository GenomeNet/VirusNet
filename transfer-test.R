# This is a test of the transfer-learning capability
# Part 1. train a model using train_type = "lm"
# Part 2. freeze layer and create a composite model for binary classificatoin
# Part 3. train a model using train_type = "label_folder"

# The LM model can be downloaded from https://genomenet.s3.eu-central-1.amazonaws.com/virus.tar.gz
# The NCBI model where the LM was trained on can be downloaded from https://genomenet.s3.eu-central-1.amazonaws.com/ncbi_genomes.tar.gz which was also used for the Bacteria set of the binary classifier 
# The viral data can be downloaded from https://genomenet.s3.eu-central-1.amazonaws.com/virus.tar.gz

library(deepG)
library(keras)

# Part 1. Run a wavenet model on SGB data
model <- create_model_lstm_cnn(
  maxlen = 400,
  layer.size = 1024,
  layers.lstm = 3,
  solver = "adam",
  learning.rate = 0.001,
  use.cudnn = TRUE,
  num_targets = 4,
  vocabulary.size = 4,
  bidirectional = TRUE)
  
summary(model)

deepG::trainNetwork(train_type = "lm", 
             model = model,
             path = "/workspace/Resources/Datasets/sgb_segata/representatives/train",
             path.val = "/workspace/Resources/Datasets/sgb_segata/representatives/validation/",
             checkpoint_path = "/workspace/Resources/Checkpoints", 
             validation.split = 0.1,
             run.name = "lstm_3",
             batch.size = 128,
             epochs = 2,
             max.queue.size = 100,
             lr.plateau.factor = 0.9,
             patience = 5,
             cooldown = 2,
             steps.per.epoch = 1000,
             step = 1000,
             vocabulary = c("A", "C", "T", "G"),
             tensorboard.log = "/workspace/Tensorboard",
             shuffleFastaEntries = TRUE,
             output = list(none = FALSE, 
                           checkpoints = TRUE, 
                           tensorboard = FALSE,
                           log = TRUE,
                           serialize_model = FALSE,
                           full_model = FALSE
             ),
             reverseComplements = FALSE,
			 wavenet_format = F)

# load the trianed model
model <- create_model_lstm_cnn(
  maxlen = 200,
  layer.size = 512,
  layers.lstm = 2,
  solver = "adam",
  learning.rate = 0.001,
  use.cudnn = TRUE,
  num_targets = 4,
  vocabulary.size = 4,
  bidirectional = FALSE)

#model <- keras::load_model_weights_hdf5(model, filepath =  '/workspace/Resources/Checkpoints/lstm_3_checkpoints/Ep.001-val_loss1.39-val_acc0.235.hdf5')

model <- keras::load_model_hdf5('/workspace/Resources/Models/ncbi_0_397.hdf5')

# remove the last 2 layers from model (dense with size 4 for vocabulary), that we can connect the internal layers with our new layers
keras::pop_layer(model)
keras::pop_layer(model)

keras::freeze_weights(model)

# create a composite model that includes the base + more layers
composite_model <- keras_model_sequential() %>%
  model %>%
  layer_dense(units = 1024, activation = "relu") %>%
  layer_dense(2) %>%
  keras::layer_activation("softmax")

optimizer <-  keras::optimizer_adam(lr = 0.01)
composite_model %>% keras::compile(loss = "categorical_crossentropy",
                           optimizer = optimizer, metrics = c("acc"))

deepG::trainNetwork(train_type = "label_folder", 
             model = composite_model,
             path = c("/workspace/Resources/Datasets/ncbi_genomes/train",
             "/workspace/Resources/Datasets/phages/all/virus"),
             path.val = c("/workspace/Resources/Datasets/ncbi_genomes/test",
             "/workspace/Resources/Datasets/phages/all/validation/virus"),
             checkpoint_path = "/workspace/Resources/Checkpoints", 
             validation.split = 0.2,
             run.name = "transfer_ncbi4",
             batch.size = 100,
             epochs = 10,
             max.queue.size = 50,
             lr.plateau.factor = 0.9,
             patience = 5,
             cooldown = 2,
             steps.per.epoch = 1000,
             step = 400,
             vocabulary = c("A", "C", "T", "G"),
             tensorboard.log = "/workspace/Tensorboard",
             shuffleFastaEntries = TRUE,
             output = list(none = FALSE, 
                           checkpoints = TRUE, 
                           tensorboard = FALSE,
                           log = TRUE,
                           serialize_model = FALSE,
                           full_model = FALSE
             ),
             labelVocabulary = c("genome", "virus"),
             reverseComplements = FALSE,
			 wavenet_format = F)
