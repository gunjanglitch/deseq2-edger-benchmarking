**DESeq2 vs edgeR Benchmark on Brown Bear RNA-seq (Hibernation Study)**
=======================================================================

This repository contains my complete analysis comparing **DESeq2** and **edgeR** on the same RNA-seq dataset to answer a commonly asked question:

> **“Why choose DESeq2 over edgeR (or vice-versa)?”**

Instead of relying on theory, I tested both tools head-to-head using real RNA-seq data.

🧬 **Dataset**
--------------
*   **Accession** [E-MTAB-14479](https://www.ebi.ac.uk/biostudies/ArrayExpress/studies/E-MTAB-14479?query=%22differential%20gene%20expression%22)

*   **Organism:** Brown Bear (_Ursus arctos_)
    
*   **Condition:** Hibernation vs Active (pre-hibernation)
    
*   **Source:** ArrayExpress raw count matrix
    
*   **Type:** Bulk RNA-seq
    
*   **Design:** Two-group comparison (Hibernation vs Awake)
    

📦 **Tools Compared**
---------------------

*   **DESeq2**
    
*   **edgeR**
    

Both were run independently from raw counts → normalization → dispersion estimation → DE analysis.

📊 **Key Results**
------------------

| Metric                                   | DESeq2 | edgeR |
|------------------------------------------|--------|-------|
| **Total DEGs (FDR < 0.05)**              | 8,964  | 8,324 |
| **Upregulated**                          | 4,463  | 4,113 |
| **Downregulated**                        | 4,501  | 4,216 |
| **Common DEGs**                          | 8,275  | —     |


### **Correlations Between Methods**

*   **log2FC correlation:** 0.9997
    
*   **Adjusted p-value correlation:** 0.9973
    
*   **Top-100 DEG overlap:** 80 / 100
    
*   **Sign concordance:** 50.86%
    
*   **rlog vs logCPM correlation:** 0.9901
    
*   **Dispersion trend correlation:** 0.93
    

Both methods showed extremely high agreement on the biological signal.

🧠 **When Does Each Tool Shine?**
---------------------------------

### ✔️ **Choose DESeq2 if your dataset has:**

*   Moderate sample size
    
*   Higher biological variation
    
*   Need for strong fold-change shrinkage
    
*   Need for slightly higher sensitivity
    
*   Preference for rlog/VST stabilization
    

### ✔️ **Choose edgeR if your dataset has:**

*   Very small sample size (n = 2–3 per group)
    
*   Many highly expressed genes
    
*   Need for conservative/strict DEG lists
    
*   Preference for speed and efficiency
    
*   Desire for a cleaner, tighter DEG set


---

## 📫 Contact

📧 Feel free to connect or reach out:  
🔗 🌐 [LinkedIn](https://www.linkedin.com/in/gunjan-sarode/) | 📫 gunjansarode.bioinfo@gmail.com  
🐛 Found a bug or have a question? Open an issue on this repo!

---

