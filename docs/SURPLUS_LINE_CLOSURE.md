# Handoff Summary — Surplus / P2 / Surrogate Line (Closed)

**Proyecto:** Cardiac CCTP / Systemic Tau + RECD  
**Línea cerrada:** Proxy redesign + surplus-over-surrogate (P2-10 + Phase α/β)  
**Fecha de cierre:** 2026-07-14  
**Estado:** Diagnostic line closed. **No** detector or I0 promotion.  
**Clinical claims:** **NONE**.

This note is the **authoritative line-close** for the local surplus research arc (diagnostic → redesign → Phase α → lag screen → Phase β). Detail lives in the linked reports; do not re-open this design without new data or a new statistic.

---

## 1. Objetivo de esta línea

Mejorar la geometría del proxy bivariado y el estadístico de surplus para reducir el piso en controles sanos (NSRDB) y recuperar una señal clara de reorganización pre-VF en la ventana de approach.

Se exploraron:

| Fase | Contenido |
|------|-----------|
| NSRDB diagnostic | Por qué \(S_t\) es ruidoso en sanos |
| Proxy redesign | Propuestas geométricas (sin tocar `build_bivariate_proxy` de producción) |
| **Phase α** | Bake-off de 5 proxies (R0, P1 residual AR, P2 lag, P3 signed, P4 residual lineal) |
| Lag screen | Lags en P2 (\(\ell \in \{5,10,20\}\)) |
| **Phase β** | Surplus-over-surrogate (phase-randomize channel 2) sobre **P2-10** |

---

## 2. Hallazgos principales (honestos)

### Logros

- **P2** (retardo \(\ell=10\)) fue el único brazo que bajó significativamente el piso en sanos (~**28%** reducción en media de \(S\) y casi eliminación de ventanas \(S \ge 0.40\)).
- Los surrogates (Phase β) recentraron el null en sanos: media de \(S^{\mathrm{ex}} \approx 0\) en NSRDB.
- El piso residual de P2-10 es en gran medida **explicable por el null de fase del canal 2** (median \(S_{\mathrm{surr}} \approx S\) en NSRDB).

### Limitaciones persistentes

- Ningún proxy ni variante de null produjo una subida sistemática en el approach respecto al basal (\(\Delta\) approach−basal sigue \(\approx 0\) o negativo).
- El contraste evento vs control se **debilitó** al aplicar surrogates (Cohen’s \(d\) bajó de ~**0.78** a ~**0.48** en \(S^{\mathrm{ex}}\)).
- El screening pre-registrado de Phase β dio **`interesting_screening = false`**.

### Conclusión de la línea

El problema del piso en sanos era principalmente **geométrico del proxy** y se resolvió en gran medida con **P2 + surrogate**.

Sin embargo, la **reorganización patológica pre-VF no se manifiesta** como aumento detectable de surplus bajo este diseño bivariado (Bandt–Pompe \(m=3\) + TV) en los datos públicos disponibles.

---

## 3. Qué se descarta / archiva

| Acción | Estado |
|--------|--------|
| Promover \(S^{\mathrm{ex}}\), \(S^{z}\) o P2 (cualquier lag) a detector o brazo I0 | **No** |
| Reabrir grids de lags, residuales P1/P4, o más nulls sobre este pipeline | **No** |
| Claims clínicos o de especificidad/FAR con estas métricas | **No** |
| Mutar `build_bivariate_proxy` de producción | **No** |
| Mezclar esta línea en el ranking detector locked (abs-z primary) | **No** |

---

## 4. Valor preservado (lo que aprendimos)

1. El acoplamiento estructural entre RR y \(|\Delta\mathrm{RR}|\) era el principal generador del piso alto en sanos.
2. Retardar el segundo canal (**P2**) es una modificación efectiva y simple para reducir ese acoplamiento simultáneo.
3. Los surrogates de fase en channel 2 son suficientes para recentrar el null en datos sanos.
4. El surplus (TV sobre patrones ordinales bivariados) tiene un **techo práctico** en estos datos públicos para detectar reorganización pre-VF.

---

## 5. Recomendación principal

**Priorizar Track B** (datos institucionales / controles no-evento device-matched).

Los controles públicos sanos (NSRDB) tienen un límite estructural para cualquier métrica de surplus de este tipo. Controles clínicos reales permitirán evaluar si la señal existe en condiciones más relevantes.

### Vías secundarias (solo si se quiere continuar explorando *otro* diseño)

- Proxies con canales más ortogonales (RR vs RMSSD local, entropía, o conteo de ectopías).
- Medidas de sinergia más cercanas al RECD Nivel 3 (excess3 / tres canales).

No reanudar micro-variantes de null o lag sobre P2-10 + TV sin un rediseño explícito.

---

## 6. Artefactos clave de esta línea

| Tipo | Path |
|------|------|
| Line close (este doc) | `docs/SURPLUS_LINE_CLOSURE.md` |
| NSRDB diagnostic | `docs/SURPLUS_NSRDB_DIAGNOSTIC.md` |
| Redesign | `docs/SURPLUS_PROXY_REDESIGN.md` |
| Phase α | `docs/SURPLUS_PROXY_BAKEOFF.md` |
| Lag screen | `docs/SURPLUS_P2_LAG_SCREEN.md` |
| Phase β | `docs/SURPLUS_P2_SURROGATE_PHASEB.md` |
| Código | `code/proxy_variants.py`, `code/surplus_surrogates.py`, runners `run_surplus_*` |
| Tests | `tests/test_proxy_variants.py`, `tests/test_surplus_surrogates.py` |
| Resultados | `results/surplus_*` (diag, bakeoff, p2lag, p2_surrogate; full cohort, no smoke) |

**Shipped production:** unchanged — **v1.4.0** @ `6f19efb`; ranking abs-z primary.  
**Local stack:** uncommitted archive candidate (Track D if user asks).

Todos los experimentos fueron **diagnósticos**, reproducibles y **sin claims clínicos**.

---

## Cierre de la línea

Esta vía de investigación (rediseño de proxy + surplus bivariado + nulls) ha sido explorada de forma exhaustiva y honesta. Los resultados son claros y se **archivan** para referencia futura.

**Successor restart:** leer `HANDOFF.md` + este archivo. Preferir Track **B** (institucional). Confirmar antes de cualquier commit. No reabrir I0 / R0–R5 / P1–P4 / continuous lags / \(S^{\mathrm{ex}}\) como detector.
