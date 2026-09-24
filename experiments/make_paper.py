import json, os, sys
sys.path.insert(0, "/home/sandbox/mega27/paperlib")
from paper import build_paper

R = json.load(open(os.path.join(os.path.dirname(__file__), "..", "results", "results.json")))
CV = json.load(open(os.path.join(os.path.dirname(__file__), "..", "results", "cv_results.json")))
fig = os.path.join(os.path.dirname(__file__), "..", "results", "figures", "gse2034.png")

build_paper(
    os.path.join(os.path.dirname(__file__), "..", "paper",
                 "MEGA27-16-gse2034-relapse-ml.docx"),
    "Honest machine learning on a real oncology cohort: bone-relapse "
    "prediction from tumor expression in GEO GSE2034",
    "Udita Phookan - MEGA-PROGRAM-27, item 16 (real-data ML study)",
    "We apply deep and classical models to a real public oncology dataset "
    f"(GSE2034: {R['n_samples']} breast-tumor expression profiles, "
    f"{R['n_genes']} probes, bone-relapse labels, prevalence "
    f"{R['prevalence']:.2f}) under leakage-free protocols. Differential "
    f"expression finds {R['n_de_fdr05']} genes at FDR<0.05. A 70/30 split "
    f"gives CNN {R['cnn_acc']:.3f} and logistic {R['logreg_acc']:.3f} "
    f"against a majority baseline of {R['majority_baseline']:.3f} - a "
    "negative we preserve rather than fish past. A 5-fold cross-validated "
    "pivot with fold-contained feature selection yields logistic AUC "
    f"{CV['logreg_auc'][0]:.3f} +/- {CV['logreg_auc'][1]:.3f} - signal "
    "above chance but below clinical utility, reported as such.",
    [
        ("Introduction", [
            "Expression signatures for relapse are a classic "
            "computational-oncology target (Wang et al. 2005 built a "
            "76-gene signature on this cohort). We ask a narrower, "
            "honestly evaluated question: with strict fold-contained "
            "feature selection and no test leakage, how much relapse "
            "signal is recoverable, and do deep models help?",
            "Protocol gates (locked before evaluation): (i) feature "
            "selection inside each training fold only; (ii) report "
            "majority-baseline alongside accuracy; (iii) preserve "
            "negative results.",
        ]),
        ("Data and methods", [
            "GSE2034 series matrix (GEO, real download): log2 transform, "
            "per-gene z-score. Labels: bone relapse (1=yes, 0=no). "
            "Differential expression: Welch t per gene, Benjamini-"
            "Hochberg FDR 0.05. Models: 1D-CNN over the expression vector "
            "(two conv layers, width 32) and from-scratch multinomial "
            "logistic regression; AUC by Mann-Whitney.",
        ]),
        ("Results - pass 1 (negative, preserved)", [
            f"70/30 split, top-400 panel: CNN {R['cnn_acc']:.3f}, logistic "
            f"{R['logreg_acc']:.3f}, majority {R['majority_baseline']:.3f}. "
            "Neither model beats the baseline: the first-pass pipeline "
            "does not support a predictive claim.",
        ]),
        ("Results - pivot (5-fold CV)", [
            f"Top-100 panel per fold, 5 folds: logistic accuracy "
            f"{CV['logreg_acc'][0]:.3f} +/- {CV['logreg_acc'][1]:.3f}, "
            f"AUC {CV['logreg_auc'][0]:.3f} +/- {CV['logreg_auc'][1]:.3f}; "
            f"CNN accuracy {CV['cnn_acc'][0]:.3f} +/- "
            f"{CV['cnn_acc'][1]:.3f}; majority {CV['majority'][0]:.3f}. "
            "AUC above 0.6 indicates real but weak signal; accuracy at or "
            "below majority underlines that point estimates on 286 samples "
            "are fragile. Verdict: signal-above-chance, not "
            "clinic-ready - consistent with the literature's difficulty "
            "on this cohort.",
        ]),
        ("Limitations", [
            "Single cohort, no external validation set; probe-level "
            "features without gene-set aggregation; CNN capacity exceeds "
            "what 286 samples can constrain. The negative first pass is "
            "reported with equal prominence.",
        ]),
    ],
    figures=[(fig, "Figure 1. Left: held-out accuracy, pass 1. Right: "
              "volcano plot, DE genes at FDR<0.05 in red.")],
    tables=[("Table 1. Performance summary.",
             ["protocol", "model", "metric", "value"],
             [["70/30", "CNN", "accuracy", f"{R['cnn_acc']:.3f}"],
              ["70/30", "logistic", "accuracy", f"{R['logreg_acc']:.3f}"],
              ["5-fold CV", "logistic", "AUC", f"{CV['logreg_auc'][0]:.3f}"],
              ["5-fold CV", "CNN", "accuracy", f"{CV['cnn_acc'][0]:.3f}"],
              ["both", "majority", "accuracy", f"{R['majority_baseline']:.3f}"]])],
    references=[
        "Wang Y. et al. Gene-expression profiles to predict distant "
        "metastasis of lymph-node-negative primary breast cancer. Lancet "
        "2005;365:671-679 (GSE2034).",
        "Benjamini Y., Hochberg Y. Controlling the false discovery rate. "
        "JRSS-B 1995;57:289-300.",
        "GEO GSE2034 series matrix: ftp.ncbi.nlm.nih.gov/geo/series/"
        "GSE2nnn/GSE2034/matrix/.",
    ])
print("paper 16 written")
