# Data README

**Author:** Anas Mourahib  
**Date:** 2025-07-08

This README describes the datasets used in **Sections 5.1 and 5.2** to [2].

---

## 📘 `data_Danube.csv`

This dataset contains river discharge values from Danube stations **7, 18, 24, 27, 30**, as used in Section 5.1.

- The original data consists of **daily discharges** between the summer months of **1960 to 2010**.
- To remove temporal dependence, the series was **declustered** following the approach in [1].

---

## 📘 `10_Industry_Portfolios_Daily.csv`

This dataset contains **value-weighted returns** for **10 industry portfolios**, used in Section 5.2.

- Source: [Kenneth French Data Library](https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/Data_Library/det_10_ind_port.html)
- Universe: Stocks listed on the **NYSE**, **AMEX**, and **NASDAQ**
- Time range: **July 1926 to March 2025**
- Sample size: **n = 51,922 daily observations**

### Portfolio definitions:
| ID | Industry        |
|----|-----------------|
| 1  | Nondurables     |
| 2  | Durables        |
| 3  | Manufacturing   |
| 4  | Energy          |
| 5  | HiTech          |
| 6  | Telecom         |
| 7  | Shops           |
| 8  | Health          |
| 9  | Utilities       |
| 10 | Other           |

As described in [2] and implemented in `Functions/portfolios`, we:
- Take **negative log-returns**
- Fit a **GARCH(1,1)** model to remove time dependence

---

## 📚 References

[1] Asadi, Peiman, Anthony C. Davison, and Sebastian Engelke.  
*"Extremes on river networks."* *The Annals of Applied Statistics*, 2015, pp. 2023–2050.

[2] Mourahib, Anas, Anna Kiriliouk, and Johan Segers.  
*A penalized least squares estimator for extreme-value mixture models.*  
arXiv preprint [arXiv:2506.15272](https://arxiv.org/abs/2506.15272), 2025.
