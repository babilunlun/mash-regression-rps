# Code supplement — MANIFEST

R 4.4.2（部分 Python 3 脚本以 python3 运行）。环境见 `sessionInfo/`（主分析期与 v2 补强期两份）。
各脚本按主图归属分文件夹；`*_v2/v3` 为迭代版，最终输出版本以各图 README 记录为准。

## figure1_fig2/
| 脚本 | 产出 |
|---|---|
| fig1_stratified.R | v2 补强：干预分层 DE（RYGB/生活方式）、交互检验、基线响应者分析（Fig1-S4；Figure1/tables/de_stratified_*、baseline_*） |
| figs_supp.R | v2 补强补充图 Fig1-S4 / Fig2-S8 / Fig6-S5 |
| upstream_regulators.R | v2 补强：emergent/reversal/down-250 的 TF 靶基因富集（ChEA2022+TRRUST2019，Fig2-S9） |
| make_fig1_v2.R | Fig 1 主图重出（v3 布局） |
| fig2s5_mash_specific.R / promote_fig2_gh.R / redo_fig2_fig4.R | Fig 2 MASH 特异层、panel g/h 提升重出、Fig2/Fig4 重出 |

## figure3/
| 脚本 | 产出 |
|---|---|
| fig3_mouse_sc.R / fig3_mouse_curate.R / fig3_mouse_score.R | 小鼠 sc（GSE305330，CellBender 过滤后）处理、注释、程序评分 |
| fig3_human_sc.R / fig3_human_curate*.R / fig3_human_label_score_rebuild.R / run_rebuild.sh | 人源 sc（GSE298719）重建、注释、评分（70,424 细胞版） |
| fig3_cellchat*.R / fig3_rctd.R / fig3_geomx.R / fig3_misty.R / fig3_cosmx_lr.R | CellChat、Visium RCTD、GeoMx、MISTy、CosMx 空间分析 |
| fig3_samac_cont.R / fig3s9_samac_leaveout_human.R | SAMac 连续分析、leave-out 循环性检验（Fig3-S8/S9） |
| fig3_gse136103.R | 独立人肝 atlas 复算 |
| decontx_human.R | v2 补强：DecontX 环境 RNA 校正 + 校正后耦合/检出率/SAMac 重算 + 供体级耦合 |
| 其余 fig3_*.R | 面板绘制与组合 |

## figure4/
| 脚本 | 产出 |
|---|---|
| stage1_downloads.sh → stage5_proxy_mr.R | GWAS 下载、区域构建、工具变量、MR 主流程 |
| run_coloc.R / susie_coloc_h3.R | coloc / SuSiE 共定位 |
| run_phewas.py / summarize_phewas.R | PheWAS |
| mr_core.R / meta_tiers.R / positive_control.R | MR 核心、证据分层 meta、阳性对照 |
| make_fig4_v3.R / finalize_fig4.R | Fig 4 主图 |

## figure5/
| 脚本 | 产出 |
|---|---|
| aggregate_ot.py / gnomad_constraint.py / gtex_tissue.py / tp53inp1_paralog.py | OpenTargets 成药性、LOEUF 约束、GTEx 组织特异性、旁系同源 |
| sc_expression_8genes.R | 8 基因细胞类型表达（Fig4f/Fig5 用） |
| loo_emergent.R / fig5s5_figure.R | 靶点-程序耦合（LOO 校正）与 Fig5-S5 |
| make_fig5.R | Fig 5 主图 |

## figure6/
| 脚本 | 产出 |
|---|---|
| fig6_prep_score.R / fig6_score_gse83452.R | RPS 构建与各队列评分 |
| fig6_analysis.R / fig6_sensitivity.R | 验证统计、敏感性分析 |
| fig6_loo.R | v2 补强：RPS 留一法（移除 TP53INP1）全管线重算（Fig6-S5） |
| make_fig6.R | Fig 6 主图 |

## figure7/
| 脚本 | 产出 |
|---|---|
| make_fig7.R | Fig 7 主图绘制（PLS 桥接、分泌组、诊断分组、BARD 对照、F3 分层面板组装） |
| reassemble_fig7.R | Fig 7 面板按 v3.29 最终顺序重排 |
| f3_sensitivity.R | F3 分层统计：AUC 复算、bootstrap CI、留一 AUC、置换检验、其他基线分层（§2.8 稳健性证据） |
| fig6_secreted_subset.R | 分泌子集评分与 1,000 次零分布检验（Fig 7E、Fig. S7） |
| fig6s6_secreted_subset_fig.R | Fig. S7 图绘制 |
| gse89632_rps_bard.csv、gse193066_f3_stratum.csv、rps160_secretome_annotation.csv | Fig 7 绘图与分泌子集分析输入表 |

## rps_calculator/
| 脚本 | 产出 |
|---|---|
| app.R、R/rps_score.R | RPS 在线计算器 Shiny 源码（部署于 https://babilunlun.shinyapps.io/rps-calculator/） |

## 说明
- 数据可用性：全部公共数据集（GSE106737/83452/135251/70779/193066/269493/294630/298719/305330/136103 等），登录号见文稿 Methods。
- Fig 1 原始 RMA/DE 流程为手工复现的标准 RMA（preprocessCore 单线程重编译），步骤与参数完整记录于 `Figure1/README.md` 方法要点；早期会话脚本未归档，本包收录其全部下游分析与绘图代码。
- 另有三个步骤以中间数据或终图形式归档：Fig 7 两个绘图输入表与分泌组注释表的生成脚本未归档（数据表见 figure7/）；Fig 7 的 PLS 桥接与分泌组面板的绘图脚本未归档（终图见 Figure7 图件）。
- 计算器运行所需的参考矩阵（reference_expr.rds，7.6 MB 二进制）与示例输入未收录；可用 figure6/ 评分脚本从 GSE193066 公共数据重建。
- v2 补强分析（2026-10-15）的输入输出对照见各图 README 的 v3.8 注记与 `FIGURE_MAP.md` v3.8 节。
