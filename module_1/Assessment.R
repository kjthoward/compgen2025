library(ggplot2)
library(cluster)
library(dplyr)
library(randomForest)
library(caret)

#Q 1.1
data <- read.table("./module_1/HistoneModeVSgeneExp.txt", header = TRUE)
model <- lm(measured_log2 ~ H3k4me3, data = data)
slope <- coef(model)["H3k4me3"]
print(slope)

#Q 1.3
plot(data$H3k27me3, data$measured_log2)
ggplot(data, aes(x = H3k27me3, y = measured_log2)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue") +
  labs(title = "Gene Expression vs H3K27me3", x = "H3K27me3 Values", y = "Gene Expression (measured_log2)")

#Q2.1

data <- read.table("./module_1/leukemiaExp.txt", header = TRUE, row.names = 1)

data_t <- t(data)

variances <- apply(data, 1, var)

top_genes <- names(sort(variances, decreasing = TRUE)[1:1000])
data_filtered <- data_t[, top_genes]

sil_widths <- numeric(10)
for (k in 2:10) {
  km <- kmeans(data_filtered, centers = k, nstart = 25)
  sil <- silhouette(km$cluster, dist(data_filtered))
  sil_widths[k] <- mean(sil[, 3])
}

optimal_k <- which.max(sil_widths)
print(paste("Optimal number of clusters:", optimal_k))

#Q2.3
km_final <- kmeans(data_filtered, centers = optimal_k, nstart = 25)

pca_result <- prcomp(data_filtered, scale. = TRUE)

pca_data <- as.data.frame(pca_result$x)
pca_data$Cluster <- as.factor(km_final$cluster)

ggplot(pca_data, aes(x = PC1, y = PC2, color = Cluster)) +
  geom_point(size = 3) +
  labs(title = "PCA of Samples", x = "PC1", y = "PC2", color = "Cluster") +
  theme_minimal()

#Q3.1

data <- read.table("./module_1/CpGMeth2Age.txt", header = TRUE)

response <- data$Age
predictors <- data[, -1]

set.seed(42)
train_index <- createDataPartition(response, p = 0.8, list = FALSE)
train_data <- predictors[train_index, ]
train_labels <- response[train_index]
test_data <- predictors[-train_index, ]
test_labels <- response[-train_index]

rf_model <- randomForest(x = train_data, y = train_labels, ntree = 500, importance = TRUE)

predictions <- predict(rf_model, newdata = test_data)

r_squared <- cor(predictions, test_labels)^2
print(paste("R-squared value:", r_squared))

#Q3.3
plot_data <- data.frame(Observed = test_labels, Predicted = predictions)

ggplot(plot_data, aes(x = Observed, y = Predicted)) +
  geom_point(color = "blue", alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Observed vs Predicted Age", x = "Observed Age", y = "Predicted Age") +
  theme_minimal()
