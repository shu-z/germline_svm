mkdir /xchip/beroukhimlab/Alex/CCLE_PCAWG/Figures_code/figure_data/scripts/fuzzy_matching/TA_launch
cd /xchip/beroukhimlab/Alex/CCLE_PCAWG/Figures_code/figure_data/scripts/fuzzy_matching/TA_launch

use UGER
use R-3.5

/xchip/beroukhimlab/Alex/PCAWG/Pan/Cluster/LaunchqsubUGER_20210111 -b /xchip/beroukhimlab/Alex/CCLE_PCAWG/Figures_code/figure_data/scripts/fuzzy_matching/20220418_TApaths.sh --launch -t 30 -m 8
