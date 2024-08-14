library(tidyverse)

# Create dataset for Tools A-E, reshaped to include precision, recall, F1 score, and concordance as separate rows for each tool
data <- data.frame(
  tool = rep(c('Canoes', 'CnMOPS', 'Codex', 'ExomeDepth', 'HaarVitCNV'), each=4),
  metric = rep(c('Precision', 'Recall', 'F1 Score', 'Concordance'), times=5),
  value = c(0.40, 0.3, 0.34, 0.3,   # Tool A
            0.78, 0.5, 0.6, 0.4,    # Tool B
            0.5, 0.2, 0.28, 0.2,    # Tool C
            0.55, 0.4, 0.46, 0.3,   # Tool D
            0.85, 0.8, 0.82, 0.7),  # Tool E
  group = rep(letters[1:5], each=4)
)

# Add empty bars for spacing
empty_bar <- 3
to_add <- data.frame(matrix(NA, empty_bar*nlevels(factor(data$group)), ncol(data)))
colnames(to_add) <- colnames(data)
to_add$group <- rep(levels(factor(data$group)), each=empty_bar)
data <- rbind(data, to_add)
data <- data %>% arrange(group)
data$id <- seq(1, nrow(data))

# Get the name and the y position of each label
label_data <- data
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) / number_of_bar
label_data$hjust <- ifelse(angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

# Prepare a data frame for base lines
base_data <- data %>%
  group_by(group) %>%
  summarize(start = min(id), end = max(id) - empty_bar) %>%
  rowwise() %>%
  mutate(title = mean(c(start, end)))

# Prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c(nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1, ]

# Make the plot
p <- ggplot(data, aes(x=as.factor(id), y=value, fill=metric)) +
  
  # Add the bars
  geom_bar(aes(x=as.factor(id), y=value, fill=metric), stat="identity", alpha=0.5) +
  
  # Add val=100/75/50/25 lines
  geom_segment(data=grid_data, aes(x=end, y=0.8, xend=start, yend=0.8), colour="grey", alpha=1, size=0.3, inherit.aes=FALSE) +
  geom_segment(data=grid_data, aes(x=end, y=0.6, xend=start, yend=0.6), colour="grey", alpha=1, size=0.3, inherit.aes=FALSE) +
  geom_segment(data=grid_data, aes(x=end, y=0.4, xend=start, yend=0.4), colour="grey", alpha=1, size=0.3, inherit.aes=FALSE) +
  geom_segment(data=grid_data, aes(x=end, y=0.2, xend=start, yend=0.2), colour="grey", alpha=1, size=0.3, inherit.aes=FALSE) +
  
  # Add labels for values
  annotate("text", x=rep(max(data$id), 4), y=c(0.2, 0.4, 0.6, 0.8), label=c("0.2", "0.4", "0.6", "0.8"), color="grey", size=3, angle=0, fontface="bold", hjust=1) +
  
  # Continue bar plotting
  geom_bar(aes(x=as.factor(id), y=value, fill=metric), stat="identity", alpha=0.5) +
  ylim(-0.2, 1) +
  
  # Theme settings
  theme_minimal() +
  theme(
    legend.position = "right",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm")
  ) +
  
  # Polar coordinate transformation
  coord_polar() +
  
  # Add tool names as labels
  geom_text(data=label_data, aes(x=id, y=value+0.1, label=tool, hjust=hjust), color="black", fontface="bold", alpha=0.6, size=2.5, angle=label_data$angle, inherit.aes=FALSE) +
  
  # Add base line information
  geom_segment(data=base_data, aes(x=start, y=-0.05, xend=end, yend=-0.05), colour="black", alpha=0.8, size=0.6, inherit.aes=FALSE) +
  geom_text(data=base_data, aes(x=title, y=-0.1, label=group), hjust=c(1,1,0,0,0), colour="black", alpha=0.8, size=4, fontface="bold", inherit.aes=FALSE)

# Display the plot
p
