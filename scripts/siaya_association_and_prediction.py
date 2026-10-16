"""Rerun of the Siaya genotype-phenotype association (binomial GLMs) and XGBoost bioassay prediction
on the combined QC-passed dataset. Functions are taken unchanged from the original notebooks
(glm.ipynb and xgboost-predict-pheno.ipynb); only data loading and reporting are new.
Run: python siaya_assoc_predict.py   (env with allel, statsmodels, scikit-learn, xgboost)"""
import os, sys, json, warnings
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
RES = "/mnt/user_shares/snagi/lstm_scratch/network_scratch/ampseq-agvampir-combined/results"
OUT = os.path.expanduser("~/siaya-rerun"); os.makedirs(OUT, exist_ok=True)
import pandas as pd
import numpy as np
import allel

def vcf_to_siaya_glm_data(vcf_path, df_samples, max_missing_filter=0.20, min_maf_filter=0.02 ,split_multiallelic=True, convert_genotypes=True):
    # Read VCF and create a dictionary
    vcf_df = vcf_to_siaya_df(vcf_path, df_samples, max_missing_filter=max_missing_filter, min_maf_filter=min_maf_filter)
    samples = vcf_df.columns[6:]

    if split_multiallelic:
        vcf_df = split_rows_with_multiple_alleles(vcf_df, samples)

    # Convert genotypes to 0/1/2
    if convert_genotypes:
        vcf_df = convert_genotype_to_alt_allele_count(vcf_df, samples)

    vcf_df = vcf_df.assign(snp_id=lambda x: "snp_" + x.CHROM + ":" + x.POS.astype(str) + "_" + x.ALT)
    vcf_df = vcf_df.set_index('snp_id')
    vcf_df = vcf_df.drop(columns=['CHROM', 'POS', 'FILTER_PASS', 'REF', 'ALT', 'ANN']).T
    vcf_df = pd.concat([df_samples.set_index('sample_id').query('location == "Siaya"'), vcf_df], axis=1)
    
    return vcf_df
    

def vcf_to_siaya_df(vcf_path, df_samples, max_missing_filter=0.20, min_maf_filter=0.02, query=None):
    # Read VCF file
    vcf_dict = allel.read_vcf(vcf_path, fields='*')
    samples = vcf_dict['samples']
    contig = vcf_dict['variants/CHROM'] 
    pos = vcf_dict['variants/POS']
    filter_pass = vcf_dict['variants/FILTER_PASS']
    ref = vcf_dict['variants/REF']
    alt = np.apply_along_axis(lambda x: ','.join(x[x != '']), 1, vcf_dict['variants/ALT'])
    ann = vcf_dict['variants/ANN']
    geno = allel.GenotypeArray(vcf_dict['calldata/GT'])
    
    print(f"Initial number of samples: {len(samples)}")
    print(f"Initial number of SNPs: {geno.shape[0]}")
    
    # filter to Siaya data
    siaya_mask = df_samples.eval("sample_id.str.contains('Siaya')", engine='python')
    sample_mask = np.isin(samples, df_samples[siaya_mask].sample_id)

    geno = geno.compress(sample_mask, axis=1)
    samples = samples[sample_mask]

    if query is not None:
        mask = df_samples.eval(query, engine='python')
        geno = geno.compress(mask, axis=1)
        samples = samples[mask]
    
    print(f"Final number of samples: {len(samples)}")
        
    # Calculate allele counts and apply filters
    ac = geno.count_alleles(max_allele=3)
    seg_mask = ac.is_segregating()
    miss_mask = (geno.is_missing().sum(axis=1) <= int(geno.shape[1]*max_missing_filter))
    maf_mask = ac.to_frequencies()[:, 1:].sum(axis=1) > min_maf_filter
    
    print(f"invariant SNPs removed: {np.sum(~seg_mask)}")
    print(f"SNPs removed due to high missingness (>{max_missing_filter*100}%): {np.sum(~miss_mask)}")
    print(f"SNPs removed due to low MAF (<{min_maf_filter}): {np.sum(~maf_mask)}")
    
    # Apply all filters
    mask = np.logical_and(seg_mask, np.logical_and(miss_mask, maf_mask))
    
    # Compress data based on mask
    geno = geno.compress(mask, axis=0)
    contig = contig[mask]
    pos = pos[mask]
    filter_pass = filter_pass[mask]
    ref = ref[mask]
    alt = alt[mask]
    ann = ann[mask]
    
    print(f"Final number of SNPs after all filters: {geno.shape[0]}")
    
    # Create DataFrame versions of VCF and genotypes
    vcf_df = pd.DataFrame({'CHROM': contig, 'POS': pos, 'FILTER_PASS': filter_pass, 'REF': ref, 'ALT': alt, 'ANN': ann})
    geno_df = pd.DataFrame(geno.to_gt().astype(str), columns=samples)
    
    # Combine VCF and genotype DataFrames
    vcf = pd.concat([vcf_df, geno_df], axis=1)
    
    print(f"Final DataFrame shape: {vcf.shape}")
    
    return vcf

def split_rows_with_multiple_alleles(df, samples):
    # Create an empty list to store the new rows
    new_rows = []
    # Iterate through each row
    for index, row in df.iterrows():
        alt_alleles = row['ALT'].split(',')
        # Check if there are multiple alleles in the ALT field
        if len(alt_alleles) > 1:
            for allele_num, allele in enumerate(alt_alleles):
                # Create a new row for each allele
                new_row = row.copy()
                new_row['ALT'] = allele
                # Update genotype fields
                for col in samples:
                    genotype = row[col]
                    # Split the genotype and process it
                    if genotype != './.':
                        gt_alleles = genotype.split('/')
                        new_gt = ['0' if (int(gt) != allele_num + 1 and gt != '0') else gt for gt in gt_alleles]
                        new_row[col] = '/'.join(new_gt)
                new_rows.append(new_row)
        else:
            new_rows.append(row)
    
    new_df = pd.DataFrame(new_rows).reset_index(drop=True)
    return new_df


def convert_genotype_to_alt_allele_count(df, samples):
    nalt_df = df.copy()
    # Iterate through each row
    for index, row in df.iterrows():
        # Update genotype fields
        for col in samples:
                genotype = row[col]
                if genotype != './.':
                    # Split the genotype and count non-zero alleles
                    alleles = genotype.split('/')
                    alt_allele_count = sum([1 for allele in alleles if allele != '0'])
                    nalt_df.at[index, col] = alt_allele_count
                else:
                    nalt_df.at[index, col] = np.nan

    return nalt_df
def calculate_pseudo_r2(results):
    """
    Calculate different pseudo R² measures for a fitted GLM model
    
    Parameters:
    results : statsmodels GLMResults object
        The fitted model results
        
    Returns:
    dict : Dictionary containing different pseudo R² measures
    """
    ll_null = results.null_deviance / -2
    ll_model = results.deviance / -2
    n = results.nobs
    
    r2_mcfadden = 1 - (ll_model / ll_null)
    r2_coxsnell = 1 - np.exp((2/n) * (ll_null - ll_model))
    r2_nagelkerke = r2_coxsnell / (1 - np.exp((2/n) * ll_null))
    
    return {
        'McFadden R²': r2_mcfadden,
        'Nagelkerke R²': r2_nagelkerke
    }

def glm_all_snps(snp_df):
    import statsmodels.formula.api as smf
    import statsmodels.api as sm
    from statsmodels.stats.multitest import fdrcorrection
    
    dfs = []
    for snp in snp_df.filter(like="snp_").columns:
        glm_data = snp_df[['phenotype', snp]].dropna()
        glm_data[snp] = glm_data[snp].astype(int)
        glm_data.columns = ['phenotype', 'snp']
    
        df_contingency = pd.crosstab(glm_data['phenotype'], glm_data['snp'])
        
        formula = f'phenotype ~ snp'
        logit_model = smf.glm(formula=formula, data=glm_data, family=sm.families.Binomial())
        log_results = logit_model.fit()
        
        # Calculate pseudo R² measures
        r2_values = calculate_pseudo_r2(log_results)
    
        pval_data = results_summary_to_dataframe(log_results)
        # Add R² values to the results
        for r2_name, r2_value in r2_values.items():
            pval_data[r2_name] = r2_value
            
        dfs.append(pval_data.assign(snp=snp))
        
    df_eff = pd.concat(dfs).query("index != 'Intercept'")
    res = fdrcorrection(df_eff['pvals'], alpha=0.05)
    df_eff['fdr'] = res[1]
    df_eff['fdr_sig'] = res[0]
    return process_effect_sizes(df_eff)


def process_effect_sizes(df_eff):
    df_eff = df_eff.assign(CHROM=lambda x: x.snp.str.split(":").str[0].str.replace("snp_", ""), 
                   POS_ALT=lambda x: x.snp.str.split(":").str[1]).assign(POS=lambda x: x.POS_ALT.str.split("_").str[0].astype(int))
    df_eff = df_eff.merge(snp_data[['CHROM', 'POS', 'ANN']])
    cols = ['snp'] + [col for col in df_eff.columns if col != 'snp']
    df_eff = df_eff[cols].drop(columns=['CHROM', 'POS'])
    df_eff = df_eff.set_index('snp')
    return df_eff    


def results_summary_to_dataframe(results):
    '''take the result of an statsmodel results table and transforms it into a dataframe'''
    pvals = results.pvalues
    coeff = results.params
    conf_lower = results.conf_int()[0]
    conf_higher = results.conf_int()[1]

    results_df = pd.DataFrame({"pvals":pvals,
                               "odds_ratio":np.exp(coeff),
                               "conf_lower":np.exp(conf_lower),
                               "conf_higher":np.exp(conf_higher)
                                })
    
    results_df.loc[:, 'sig'] = [True if pval <= 0.05 else False for pval in pvals]
    results_df = results_df[["odds_ratio","pvals","conf_lower","conf_higher", "sig"]]
    return results_df
import pandas as pd
import numpy as np

def compute_snp_correlations(df, correlation_threshold=0.5):
    corr_matrix = df.corr()
    mask = np.triu(np.ones_like(corr_matrix, dtype=bool))
   
    high_corr_pairs = []
    for i in range(len(corr_matrix.columns)):
        for j in range(i+1, len(corr_matrix.columns)):
            if abs(corr_matrix.iloc[i, j]) > correlation_threshold:
                high_corr_pairs.append((corr_matrix.columns[i], corr_matrix.columns[j], corr_matrix.iloc[i, j]))
    
    print(f"Total SNPs analyzed: {len(df.columns)}")
    print(f"Number of highly correlated SNP pairs (|r| > {correlation_threshold}): {len(high_corr_pairs)}")
    
    return corr_matrix, high_corr_pairs

def _unused_plot(corr_matrix):
    plt.figure(figsize=(12, 10))
    sns.heatmap(corr_matrix, annot=False, cmap='coolwarm', vmin=-1, vmax=1, center=0)
    plt.title('SNP Correlation Heatmap')
    plt.show()

def select_representative_snps(df, high_corr_pairs):
    to_remove = set()
    pos_corr_removed = 0
    neg_corr_removed = 0
    
    for snp1, snp2, corr in high_corr_pairs:
        if snp1 not in to_remove and snp2 not in to_remove:
            if corr < 0:  # Negative correlation
                if df.loc[snp1, 'odds_ratio'] < 1:
                    to_remove.add(snp1)
                    neg_corr_removed += 1
                elif df.loc[snp2, 'odds_ratio'] < 1:
                    to_remove.add(snp2)
                    neg_corr_removed += 1
                else:
                    if df.loc[snp1, 'fdr'] < df.loc[snp2, 'fdr']:
                        to_remove.add(snp2)
                        neg_corr_removed += 1
                    else:
                        to_remove.add(snp1)
                        neg_corr_removed += 1
            else:  # Positive correlation
                if df.loc[snp1, 'fdr'] < df.loc[snp2, 'fdr']:
                    to_remove.add(snp2)
                    pos_corr_removed += 1
                else:
                    to_remove.add(snp1)
                    pos_corr_removed += 1
    
    selected_snps = np.array(list(set(df.index) - to_remove))
    
    print(f"Total SNPs before selection: {len(df)}")
    print(f"SNPs removed due to positive correlation: {pos_corr_removed}")
    print(f"SNPs removed due to negative correlation: {neg_corr_removed}")
    print(f"Total SNPs removed: {len(to_remove)}")
    print(f"SNPs retained after selection: {len(selected_snps)}")
    
    return selected_snps

def linked(df_genos, df_eff, correlation_threshold=0.5):
    corr_matrix, high_corr_pairs = compute_snp_correlations(df_genos.filter(like='snp_'), correlation_threshold)
    selected_snps = select_representative_snps(df_eff, high_corr_pairs)
    return np.sort(selected_snps)
from sklearn.model_selection import cross_val_score, cross_val_predict, GridSearchCV, StratifiedKFold
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.dummy import DummyClassifier
import xgboost as xgb
from scipy.stats import uniform, randint

def prepare_data(df):
    """
    Prepare the data for modeling by:
    1. Dropping unnecessary columns
    2. Encoding categorical variables
    3. Handling missing values
    """
    # Create copy of dataframe
    df = df.copy()
    
    # Drop non-feature columns
    columns_to_drop = ['taxon', 'aim_taxon', 'pca_taxon', 'tree_taxon',  'sample_name', 'index', 'index2', 'location', 'country', 'plate', 'well_letter', 
                      'well_number', 'latitude', 'longitude']
    df = df.drop([col for col in columns_to_drop if col in df.columns], axis=1)
    
    # Encode phenotype
    le = LabelEncoder()
    df['phenotype'] = le.fit_transform(df['phenotype'])
    
    # Fill missing values with mode of column
    df = df.fillna(df.mode().iloc[0])
    
    return df, le

def get_random_baseline(X, y, cv):
    """
    Calculate baseline accuracy using different random strategies with cross-validation
    """
    print("\nRandom Baseline Performance:")
    print("--------------------------")
    
    strategies = ['stratified', 'most_frequent', 'uniform']
    results = {}
    
    for strategy in strategies:
        # Initialize dummy classifier
        dummy = DummyClassifier(strategy=strategy, random_state=42)
        
        # Get cross-validation scores
        scores = cross_val_score(dummy, X, y, cv=cv, scoring='accuracy')
        
        # Store and print results
        results[strategy] = {
            'mean': scores.mean(),
            'std': scores.std()
        }
        print(f"{strategy.capitalize()} random baseline accuracy: {scores.mean():.3f} (+/- {scores.std()*2:.3f})")
    
    return results

from sklearn.model_selection import GridSearchCV
import xgboost as xgb

def optimize_xgboost_hyperparameters(X, y, cv):
    """
    Perform hyperparameter optimization using GridSearchCV
    """
    # param_grid = {
    #     'n_estimators': [100, 500],
    #     'max_depth': [3, 7, 12],
    #     'learning_rate': [0.01, 0.1, 0.3],
    #     'subsample': [0.7, 1.0],
    #     'colsample_bytree': [0.7, 1.0],
    #     'min_child_weight': [1, 5],
    #     'gamma': [0.0, 0.3],
    #     'reg_alpha': [0.0, 1.0],
    #     'reg_lambda': [1.0, 4.0]
    # }

    param_grid = {
        'n_estimators': [100, 500],
        'max_depth': [3, 7, 12],
        'learning_rate': [0.01, 0.1, 0.3],
        'subsample': [0.7, 1.0],
        'colsample_bytree': [0.7, 1.0],
        'min_child_weight': [1, 5],
        'gamma': [0.0, 0.3],
        'reg_alpha': [0.0, 1.0],
        'reg_lambda': [1.0, 4.0],
    }
    
    # Initialize base model
    base_model = xgb.XGBClassifier(
        objective='binary:logistic',
        random_state=42, n_jobs=1,
    )
    
    # Initialize GridSearchCV
    grid_search = GridSearchCV(
        estimator=base_model,
        param_grid=param_grid,
        scoring='accuracy',
        cv=cv,
        n_jobs=-1,
        verbose=1
    )
    
    # Fit GridSearchCV
    grid_search.fit(X, y)
    
    # Print results
    print("\nHyperparameter Optimization Results:")
    print("----------------------------------")
    print(f"Best cross-validation accuracy: {grid_search.best_score_:.3f}")
    print("\nBest parameters:")
    for param, value in grid_search.best_params_.items():
        print(f"{param}: {value}")
        
    return grid_search.best_estimator_


def evaluate_model_cv(model, X, y, cv, label_encoder=None):
    """
    Evaluate model using cross-validation
    """
    # Get cross-validation scores
    cv_scores = cross_val_score(model, X, y, cv=cv, scoring='accuracy')
    
    # Get predictions for detailed metrics
    y_pred = cross_val_predict(model, X, y, cv=cv)
    
    print("\nCross-Validation Results:")
    print("------------------------")
    print(f"Accuracy: {cv_scores.mean():.3f} (+/- {cv_scores.std()*2:.3f})")
    
    print("\nClassification Report:")
    print("---------------------")
    if label_encoder:
        target_names = label_encoder.classes_
        print(classification_report(y, y_pred, target_names=target_names))
    else:
        print(classification_report(y, y_pred))
    
    print("\nConfusion Matrix:")
    print("----------------")
    print(confusion_matrix(y, y_pred))
    
    return cv_scores

def get_feature_importance(model, X, y, feature_names):
    """
    Get and print feature importance from the model
    """
    # Fit model on full dataset to get feature importance
    model.fit(X, y)
    importance = model.feature_importances_
    indices = np.argsort(importance)[::-1]

    print("\nTop 10 Most Important Features:")
    print("------------------------------")
    for f in range(min(10, len(feature_names))):
        print("%d. %s (%f)" % (f + 1, feature_names[indices[f]], importance[indices[f]]))
        
    return pd.DataFrame({
        'Feature': feature_names,
        'Importance': importance
    }).sort_values('Importance', ascending=False)

def predict_phenotype(df, n_folds=5, n_iter=100):
    """
    Main function to run the entire prediction pipeline with cross-validation
    Returns all necessary variables for visualization
    """
    # Prepare data
    df, label_encoder = prepare_data(df)

    # Split features and target
    X = df.drop('phenotype', axis=1)
    y = df['phenotype']
    
    # Create stratified k-fold cross-validator
    cv = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=42)

    # Optimize hyperparameters
    print("\nOptimizing hyperparameters... This may take a while...")
    best_model = optimize_xgboost_hyperparameters(X, y, cv=cv)
        
    # Return all variables needed for visualization
    results = {
        'model': best_model,
        'X': X,
        'y': y,
        'cv_scores': evaluate_model_cv(best_model, X, y, cv, label_encoder),
        'importance_df': get_feature_importance(best_model, X, y, X.columns),
        'label_encoder': label_encoder,
        'baseline_results': get_random_baseline(X, y, cv),
        'cv': cv
    }
    
    return results

# ---------------- data ----------------
meta = pd.read_csv(f"{RES}/config/metadata.qcpass.tsv", sep="\t")
meta["phenotype"] = meta.sample_id.str.extract(r"_(Dead|Alive)_")[0]
siaya = meta[meta.sample_id.str.startswith("Siaya_Delta") & meta.phenotype.notna()].copy()
siaya["location"] = "Siaya"
print("phenotyped QC-passed Siaya samples:", len(siaya), siaya.phenotype.value_counts().to_dict())
plates = pd.crosstab(siaya.plate, siaya.phenotype); print(plates); plates.to_csv(f"{OUT}/plate_by_phenotype.csv")

snp_data = pd.read_excel(f"{RES}/vcfs/amplicons/ag-vampir-combined-snps.xlsx")
snp_data = snp_data[["CHROM","POS","ANN"]].drop_duplicates(["CHROM","POS"])  # combined table lists some positions more than once
print("snp table rows:", len(snp_data))
df_genos = vcf_to_siaya_glm_data(f"{RES}/vcfs/amplicons/ag-vampir-combined.annot.vcf", siaya, max_missing_filter=0.20, min_maf_filter=0.02)
# GLM response: survival (Alive = 1), so OR > 1 means the alternate allele is associated with survival
df_genos["phenotype"] = (df_genos["phenotype"] == "Alive").astype(int)
snps = [c for c in df_genos.columns if c.startswith("snp_")]
print("samples:", df_genos.shape[0], "SNPs after filters:", len(snps))
df_eff = glm_all_snps(df_genos)
df_eff.to_csv(f"{OUT}/glm_effects.csv")
sig = df_eff[df_eff.fdr_sig].copy()
print("tested:", len(df_eff), "FDR-significant:", len(sig))
print(sig[["odds_ratio","conf_lower","conf_higher","fdr"]].round(3).to_string())
sel = linked(df_genos[["phenotype"] + snps].dropna().astype(int), df_eff)
print("FDR-significant after LD pruning (representatives):", len(set(sig.index) & set(sel)))

# ---------------- XGBoost (original procedure) ----------------
df_x = df_genos[["phenotype"] + snps].copy()
df_x["phenotype"] = np.where(df_x.phenotype == 1, "Alive", "Dead")
res = predict_phenotype(df_x, n_folds=5, n_iter=10)
best = res["model"]; X, y, cv, le = res["X"], res["y"], res["cv"], res["label_encoder"]
from sklearn.metrics import classification_report, f1_score, accuracy_score
yp = cross_val_predict(best, X, y, cv=cv)
cr = classification_report(y, yp, target_names=le.classes_, output_dict=True)
print("\nORIGINAL PROCEDURE (grid search and CV on the same folds):")
print("grid best CV accuracy: %.3f" % res["cv_scores"].mean(), "+/-", round(res["cv_scores"].std(), 3))
print("accuracy %.3f  weighted F1 %.3f" % (accuracy_score(y, yp), f1_score(y, yp, average="weighted")))
print(pd.DataFrame(cr).T.round(2).to_string())
print("best params:", best.get_params())
print(res["importance_df"].head(10).to_string())
res["importance_df"].to_csv(f"{OUT}/xgb_importance.csv")

# ---------------- nested CV (no leakage from tuning) ----------------
from sklearn.model_selection import StratifiedKFold, GridSearchCV
outer = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
preds = np.zeros(len(y), dtype=int); accs = []
for tr, te in outer.split(X, y):
    gs = GridSearchCV(xgb.XGBClassifier(objective="binary:logistic", random_state=42, n_jobs=1),
        {'n_estimators':[100,500],'max_depth':[3,7,12],'learning_rate':[0.01,0.1,0.3],'subsample':[0.7,1.0],
         'colsample_bytree':[0.7,1.0],'min_child_weight':[1,5],'gamma':[0.0,0.3],'reg_alpha':[0.0,1.0],'reg_lambda':[1.0,4.0]},
        scoring="accuracy", cv=StratifiedKFold(3, shuffle=True, random_state=1), n_jobs=-1)
    gs.fit(X.iloc[tr], y.iloc[tr]); p = gs.predict(X.iloc[te]); preds[te] = p; accs.append(accuracy_score(y.iloc[te], p))
print("\nNESTED CV: accuracy %.3f (+/- %.3f)  weighted F1 %.3f" % (np.mean(accs), np.std(accs), f1_score(y, preds, average="weighted")))
print(pd.DataFrame(classification_report(y, preds, target_names=le.classes_, output_dict=True)).T.round(2).to_string())
json.dump({"n": int(len(y)), "orig_cv_acc": float(res["cv_scores"].mean()), "nested_acc": float(np.mean(accs)),
           "nested_f1": float(f1_score(y, preds, average="weighted"))}, open(f"{OUT}/summary.json", "w"))
