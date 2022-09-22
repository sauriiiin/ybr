source('/home/sbp29/R/Projects/ybr/scripts/crispr_screen/initialize.R')

temp <- dbGetQuery(conn, 'select a.*, b.orf_name, c.strain_id
              from YBR_CRISPR384_pos2coor a, YBR_CRISPR384_pos2orf_name b, YBR_CRISPR384_pos2strainid c
              where a.pos = b.pos and b.pos = c.pos')
temp$strain_id <- as.numeric(str_trunc(as.character(temp$strain_id), 3, side = 'left', ellipsis = ''))
temp$label <- paste(temp$orf_name, '\n', temp$strain_id, sep = '')

plate384 <- temp %>%
  filter(density == 384, plate_no == 1) %>%
  ggplot(aes(x = plate_col, y = plate_row)) +
  geom_tile(aes(fill = orf_name), col = 'black') +
  geom_text(aes(label = label), size = 1.25, angle = 45) +
  scale_y_reverse(breaks = seq(16,1,-1)) +
  scale_x_continuous(breaks = seq(1,24,1)) +
  scale_fill_discrete(name = '') +
  # facet_wrap(.~plate_no) +
  theme_linedraw()

ggsave(sprintf("%s/YBR_CRISPR_PLATE384.jpg",fig_path), plate384,
       height = 160, width = 280, units = 'mm',
       bg = 'white',
       dpi = 600)
