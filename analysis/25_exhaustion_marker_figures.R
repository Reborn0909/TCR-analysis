



suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr); library(readr)})
setwd(Sys.getenv("DATA_ROOT", "."))
OUT <- "figures/exhaustion"
group_colors <- c("E-DP"="#C8956B", "E-Trp"="#B05A5A", "H-Trp"="#6B8E9B")
sample_colors <- c("E-CD8DP-1"="#C8956B","E-CD8DP-2"="#D9AE8A","E-CD8Trp-1"="#B05A5A",
                   "E-CD8Trp-2"="#C27A7A","H-CD8Trp-1"="#6B8E9B","H-CD8Trp-2"="#8BA6B2")
RAMP <- c("#D4CFC5","#B5C1B4","#8C5A5A")
SAMP <- c("H-CD8Trp-1","E-CD8Trp-1","E-CD8DP-1","H-CD8Trp-2","E-CD8Trp-2","E-CD8DP-2")
theme_pub <- theme_minimal(base_size=7) +
  theme(plot.title=element_text(size=8.5,face="bold",hjust=.5,margin=margin(b=3)),
        plot.subtitle=element_text(size=6.8,hjust=.5),
        axis.title=element_text(size=7), axis.text=element_text(size=6,color="black"),
        axis.line=element_line(linewidth=.3,color="black"),
        axis.ticks=element_line(linewidth=.3,color="black"),
        legend.title=element_text(size=7), legend.text=element_text(size=6),
        legend.key.size=unit(.35,"cm"),
        panel.grid.major=element_line(linewidth=.25,color="grey90"), panel.grid.minor=element_blank(),
        panel.border=element_rect(linewidth=.3,fill=NA,color="black"))
theme_set(theme_pub)
save_panel <- function(p,name,w=5.0,h=4.4){ cairo_pdf(file.path(OUT,paste0(name,".pdf")),width=w,height=h); print(p); dev.off()
  ggsave(file.path(OUT,paste0(name,".png")),p,width=w,height=h,dpi=600); cat("  :",name,"\n") }

agg <- read_csv("results/exhaustion_marker_panel_donor_delta.csv", show_col_types=FALSE)
long <- read_csv("results/exhaustion_marker_panel_donor_expression.csv", show_col_types=FALSE)


dd <- agg %>% select(gene, category, d1, d2, mean_delta, up_both, down_both) %>%
  pivot_longer(c(d1,d2), names_to="donor", values_to="delta") %>%
  mutate(donor = recode(donor, d1="Donor 1", d2="Donor 2"),
         gene = factor(gene, levels = agg$gene[order(agg$mean_delta)]))
p <- ggplot(dd, aes(delta, gene)) +
  geom_vline(xintercept=0, linewidth=.35, colour="grey40") +
  geom_point(aes(colour=donor), size=1.9) +
  geom_segment(data=agg %>% mutate(gene=factor(gene, levels=levels(dd$gene))),
               aes(x=d1, xend=d2, y=gene, yend=gene), colour="grey55", linewidth=.35,
               inherit.aes=FALSE) +
  scale_colour_manual(values=c("Donor 1"="#8C5A5A","Donor 2"="#A87575"), name=NULL) +
  facet_grid(category ~ ., scales="free_y", space="free_y", switch="y") +
  labs(title="Exhaustion markers: CD8 Trp vs SS",
       subtitle="Donor-level mean log1p expression difference; bars join the two donors",
       x="Difference in mean log1p expression", y=NULL) +
  theme(strip.text.y.left=element_text(size=6, angle=0), strip.placement="outside",
        strip.background=element_rect(fill="grey95", colour="black", linewidth=.3),
        legend.position="top", panel.spacing.y=unit(2,"pt"))
save_panel(p, "E1_exhaustion_marker_delta", 4.8, 6.4)


hm <- long %>%
  mutate(col = factor(paste0(recode(group, "H-Trp"="SS","E-Trp"="Pre-sel","E-DP"="Post-sel"),
                            "\n", donor),
                      levels = c("SS\nDonor 1","Pre-sel\nDonor 1","Post-sel\nDonor 1",
                                 "SS\nDonor 2","Pre-sel\nDonor 2","Post-sel\nDonor 2")),
         gene = factor(gene, levels=agg$gene[order(agg$mean_delta)]))

hm <- hm %>% group_by(gene) %>%
  mutate(z = as.numeric(scale(mean_log1p))) %>% ungroup() %>%
  mutate(z = ifelse(is.finite(z), z, 0))
p <- ggplot(hm, aes(col, gene, fill=z)) +
  geom_tile(colour="white", linewidth=.35) +
  scale_fill_gradientn(colours=RAMP, name="Row z-score",
                       limits=c(-2.2, 2.2), oob=scales::squish) +
  labs(title="Exhaustion-related marker panel",
       subtitle="Mean log1p normalised expression, each gene standardised across the six libraries",
       x=NULL, y=NULL) +
  theme(axis.text.x=element_text(angle=45,hjust=1), panel.grid=element_blank())
save_panel(p, "E2_exhaustion_marker_heatmap", 4.0, 6.2)


top14 <- agg %>% filter(up_both) %>% arrange(desc(mean_delta)) %>% slice_head(n=14)
up <- top14 %>% select(gene, category, d1, d2) %>%
  pivot_longer(c(d1,d2), names_to="donor", values_to="delta") %>%
  mutate(donor = recode(donor, d1="Donor 1", d2="Donor 2"),
         gene   = factor(gene, levels = rev(top14$gene)))
p <- ggplot(up, aes(delta, gene, fill=donor)) +
  geom_col(position=position_dodge(width=.72), width=.66, colour="grey30", linewidth=.2) +
  geom_vline(xintercept=0, linewidth=.35, colour="grey40") +
  scale_fill_manual(values=c("Donor 1"="#8C5A5A","Donor 2"="#A87575"), name=NULL) +
  labs(title="Markers up in CD8 Trp in both donors",
       subtitle="Top 14 by mean difference (Pre-sel minus SS)", x="Difference in mean log1p expression", y=NULL) +
  theme(legend.position="top")
save_panel(p, "E3_upregulated_markers", 4.6, 3.4)
cat("\n\n")
