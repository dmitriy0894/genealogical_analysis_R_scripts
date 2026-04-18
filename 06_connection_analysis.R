library(openxlsx)
library(dplyr)
library(ggplot2)

# --- 1. Workspace Setup & File Selection ---
cat("Please select the IBD results file (Excel)...\n")
ibd_results.fn <- file.choose()

work_dir <- dirname(ibd_results.fn)
setwd(work_dir)
cat("Working directory set to:", work_dir, "\n")

# --- 2. Load & Prepare Data ---
ibd_data <- read.xlsx(ibd_results.fn)

# Define relationship colors
rel_colors <- c(
  "Monozygotic twin / Clone" = "purple",
  "Parent-offspring" = "red",
  "Full sib"         = "blue",
  "2nd Degree"       = "green",
  "3rd Degree"       = "orange"
)

# Define factor levels order
rel_levels <- names(rel_colors)

cat("Processing data and calculating connections...\n")

# Keep only significant connections
valid_network_links <- ibd_data %>%
  filter(kinship > 0, Relationship != "Unrelated", !is.na(Relationship))

# Prepare data for counting (ID1 and ID2)
part1 <- valid_network_links %>% select(Cultivar = ID1, Relationship)
part2 <- valid_network_links %>% select(Cultivar = ID2, Relationship)
all_nodes <- bind_rows(part1, part2)

# Ensure correct factor order for relationships
all_nodes <- all_nodes %>%
  mutate(Relationship = factor(Relationship, levels = rel_levels))

# Count connections by relationship type
rel_counts <- all_nodes %>%
  group_by(Cultivar, Relationship) %>%
  summarise(Count = n(), .groups = 'drop')

# Calculate total connections per cultivar
total_counts <- rel_counts %>%
  group_by(Cultivar) %>%
  summarise(Total = sum(Count), .groups = 'drop') %>%
  arrange(Total)

# Bind factors for sorting on the plot
rel_counts$Cultivar <- factor(rel_counts$Cultivar, levels = total_counts$Cultivar)
total_counts$Cultivar <- factor(total_counts$Cultivar, levels = total_counts$Cultivar)

# --- 3. Layout Parameters ---
max_char <- max(nchar(as.character(total_counts$Cultivar)))
left_box_width <- max_char * 0.045 + 0.35 

# --- 4. Plot Generation ---
cat("Building the plot...\n")

bar_plot <- ggplot() +
  
  # 1. Left panel (Light gray box for text)
  geom_col(data = total_counts, aes(x = -left_box_width, y = Cultivar), 
           fill = "gray95", color = "black", width = 0.65, linewidth = 0.5) +
  
  # 2. Text inside the left panel
  geom_text(data = total_counts, 
            aes(x = -left_box_width + 0.03, y = Cultivar, 
                label = paste0(Cultivar, " (edges: ", Total, ")")), 
            hjust = 0, fontface = "plain", size = 3.8, color = "black") +
  
  # 3. Colored blocks (100% stacked)
  geom_col(data = rel_counts, aes(x = Count, y = Cultivar, fill = Relationship), 
           position = position_fill(reverse = TRUE), width = 0.65, color = "black", linewidth = 0.5) +
  
  # 4. Numbers inside colored blocks
  geom_text(data = subset(rel_counts, Count > 0), 
            aes(x = Count, y = Cultivar, label = Count, group = Relationship), 
            position = position_fill(vjust = 0.5, reverse = TRUE), 
            color = "white", fontface = "plain", size = 4) +
  
  # Configure colors and legend
  scale_fill_manual(values = rel_colors) +
  labs(fill = "Edges of:") +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "plain", size = 11),
    legend.text = element_text(size = 10, face = "plain", margin = margin(l = 10, r = 20)),
    legend.spacing.x = unit(0.3, "cm"),
    plot.margin = margin(10, 5, 10, 5) 
  )

# --- 5. Export Plot (PNG, 600 DPI) ---
output_bar_png <- "Cultivars_Connections_Planks.png"
dynamic_height <- max(5, nrow(total_counts) * 0.4)

ggsave(filename = output_bar_png, plot = bar_plot, device = "png", 
       width = 9, height = dynamic_height, units = "in", dpi = 600, bg = "white")

cat("SUCCESS! Plot saved as:", output_bar_png, "\n")
