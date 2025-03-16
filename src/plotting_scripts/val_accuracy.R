library(ggplot2)

data_file_name <- "data/val_accuracy.csv"

data <- read.csv(data_file_name)

data$color <- ifelse(data$epoch == 6, "#d200ef", "black")

plot <- ggplot(data, aes(x=epoch, y=ruby.sea.3...val_acc)) +
        geom_point() + 
        geom_line(color="black", size=0.4) +
        geom_point(aes(color = color), size=0.6,  stroke = 0.9) + 
        scale_x_continuous(limits=c(0,9), breaks=seq(0, 9, by=1)) +
        scale_y_continuous(limits=c(0.9, 1)) +
        scale_color_identity() + 
        labs(title="", x="Epoch", y="Validation accuracy") +
        theme_minimal() + 
        theme(panel.grid = element_blank()) + 
        theme(axis.line = element_line(color = "black"))
plot

ggsave("plots/validation_accuracy.pdf", plot=plot, width=8, height=6)