import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.linear_model import Ridge
from sklearn.model_selection import cross_val_score, KFold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score
import matplotlib.pyplot as plt

df = pd.read_csv('data/selectivity_with_descriptors.csv')

features = ['MW', 'LogP', 'HBD', 'HBA', 'TPSA', 'RotBonds', 'ArRings', 'RingCount']
X = df[features].values
y = df['selectivity_ratio'].values

print(f"Dataset: {len(df)} compounds")
print(f"Features: {features}")

scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

cv = KFold(n_splits=5, shuffle=True, random_state=42)

models = {
    'Ridge regression': Ridge(alpha=1.0),
    'Random Forest':    RandomForestRegressor(n_estimators=100, random_state=42),
    'Gradient Boosting': GradientBoostingRegressor(n_estimators=100, random_state=42),
}

print("\n── 5-fold cross-validation R² ──")
results = {}
for name, model in models.items():
    X_in = X_scaled if name == 'Ridge regression' else X
    scores = cross_val_score(model, X_in, y, cv=cv, scoring='r2')
    results[name] = scores
    print(f"{name:25s}: R² = {scores.mean():.3f} ± {scores.std():.3f}")

# Feature importance from Random Forest
rf = RandomForestRegressor(n_estimators=100, random_state=42)
rf.fit(X, y)
importances = pd.Series(rf.feature_importances_, index=features).sort_values(ascending=False)
print("\n── Random Forest feature importances ──")
print(importances.round(3))

# Plot
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Feature importance
importances.plot(kind='bar', ax=axes[0], color='steelblue', alpha=0.8)
axes[0].set_title('Feature importance for 5-HT2A/2B selectivity')
axes[0].set_ylabel('Importance')
axes[0].tick_params(axis='x', rotation=45)

# Predicted vs actual
rf_pred = cross_val_score(rf, X, y, cv=cv, scoring='r2')
from sklearn.model_selection import cross_val_predict
y_pred = cross_val_predict(rf, X, y, cv=cv)
axes[1].scatter(y, y_pred, alpha=0.3, s=10, color='steelblue')
axes[1].plot([-3, 3], [-3, 3], 'r--', alpha=0.5)
axes[1].set_xlabel('Actual selectivity (pKi_2A - pKi_2B)')
axes[1].set_ylabel('Predicted selectivity')
axes[1].set_title(f'Random Forest: R² = {rf_pred.mean():.3f}')

plt.tight_layout()
plt.savefig('analysis/qsar_model.png', dpi=150, bbox_inches='tight')
plt.close()

print("\nSaved analysis/qsar_model.png")
