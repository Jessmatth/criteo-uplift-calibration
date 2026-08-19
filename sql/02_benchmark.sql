-- Randomized ground truth, computed on all 13,979,592 rows.
--
-- Treatment was assigned at random, so no adjustment is required or wanted: the
-- difference in outcome rates between arms IS the causal effect of the campaign.
-- Everything downstream in this repository is graded against these numbers.
--
-- Sufficient statistics only. The delta-method covariance for the instrumental
-- variable estimate needs E[y*d] per arm, which is why ydbar is carried here.
SELECT
    treatment,
    count(*)                          AS n,
    sum(conversion)                   AS conversions,
    avg(conversion)                   AS ybar_conversion,
    sum(visit)                        AS visits,
    avg(visit)                        AS ybar_visit,
    avg(exposure)                     AS dbar_exposure,
    avg(conversion * exposure)        AS ydbar_conversion,
    avg(visit * exposure)             AS ydbar_visit
FROM criteo
GROUP BY treatment
ORDER BY treatment;
