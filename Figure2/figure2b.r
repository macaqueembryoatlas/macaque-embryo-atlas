library(ggalluvial)
library(ggplot2)
library(dplyr)
seu.rna<-readRDS('/Users/mac_tian/embryo/telecephalon/data/whole(npc.neuron).rna.meta.Rds')
seu.rna$celltype_f2<-seu.rna$celltype_f

color.neuron2 <- c('Gaba_c2'='#c27e58',
                   'Undefined'="#4693A8",
                   'Gaba_CB_c1'="#c153b4",
                   'Gaba_MB_c1'='#6CA9CF',
                   'Glu_MB_c1'='#6484f5',
                   'Glu_CB_c1'='#e213a7',
                   'Gaba_MGE'='#5aae6b',
                   'Glu_CB_c2'='#dc8ac5',
                   'Glu_c1'="#F48829",
                   'Gaba_CGE'="#00800c",
                   'Gaba_MB_c2'='#256595', 
                   'Gaba_BS_c1'="#75707c", 
                   'Gaba_Purk'="#D6CA99", 
                   'Glu_MB_c2'='#3758cc',
                   'Gaba_c1'="#B15928",
                   'Glu_Tel_Newborn_latter'='#67a26d',
                   'Glu_Tel_Mature'='#009900',
                   'Glu_Tel_Newborn_earlier'="#cce0ce",
                   'Gaba_MB_c3'='#0069ff')
color.npc<-c('RG_c1'="#646f9b",
             'CB_RG'="#a2ac86",
             'MB_RG'="#d37959",
             'Prec1'="#18d600",
             'Dien_RG'="#ffbfbe",
             'MO_RG'="#a700d6",
             'Tel_RG'='#ff6311',
             'Prec2'="#ffa600",
             'IPC'="#dc9f46")
color2<-c(color.neuron2,color.npc)

counts<-seu.rna%>%group_by(Day,celltype_f2)%>%
  summarise(count=n())%>%
  ungroup()%>%
  group_by(Day)%>%
  mutate(propotion=count/sum(count))
# 修正时间点转换（假设原始格式为"DayX"而非"D"开头）
counts <- counts %>%
  mutate(Day_numeric = as.numeric(gsub("D", "", Day)))  # 正确提取数字部分

levels.celltype<-c('Tel_RG','CB_RG','Dien_RG','MB_RG','MO_RG','RG_c1','IPC','Prec1','Prec2',
                   'Glu_Tel_Newborn_earlier','Glu_Tel_Newborn_latter','Glu_Tel_Mature',
                   'Glu_MB_c1','Glu_MB_c2','Glu_CB_c1','Glu_CB_c2','Glu_c1',
                   'Gaba_MGE','Gaba_CGE','Gaba_MB_c1','Gaba_MB_c2','Gaba_MB_c3','Gaba_CB_c1','Gaba_BS_c1','Gaba_Purk','Gaba_c1','Gaba_c2')
counts$celltype_f2<-factor(counts$celltype_f2,levels=levels.celltype)



pdf('./fig2b.pdf',width=8,height=8)
ggplot(counts, aes(x = Day_numeric, y = propotion, 
                   alluvium = celltype_f2, stratum = celltype_f2)) +
  geom_flow(aes(fill = celltype_f2), alpha = 1, width = 0) +  # 使用geom_flow而非geom_alluvium
  geom_stratum(aes(fill = celltype_f2), width = 0, color = "grey") +
  scale_fill_manual(values = color2) +
  scale_x_continuous(
    breaks = unique(counts$Day_numeric),
    labels = unique(counts$Day)
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  ) +
  labs(title = "scRNA-seq data",
       x = "Day", y = "Proportion", fill = "Celltype")+
  coord_fixed(ratio=20)
dev.off()
