-- Stage the raw archive into DuckDB.
--
-- read_csv_auto reads the gzip stream directly: the 3.5 GB uncompressed CSV is
-- never written to disk. Columns are typed explicitly rather than sniffed so
-- that a change in the upstream file surfaces as a load error instead of a
-- silently different schema.
CREATE OR REPLACE TABLE criteo AS
SELECT
    f0::DOUBLE   AS f0,
    f1::DOUBLE   AS f1,
    f2::DOUBLE   AS f2,
    f3::DOUBLE   AS f3,
    f4::DOUBLE   AS f4,
    f5::DOUBLE   AS f5,
    f6::DOUBLE   AS f6,
    f7::DOUBLE   AS f7,
    f8::DOUBLE   AS f8,
    f9::DOUBLE   AS f9,
    f10::DOUBLE  AS f10,
    f11::DOUBLE  AS f11,
    treatment::TINYINT  AS treatment,
    conversion::TINYINT AS conversion,
    visit::TINYINT      AS visit,
    exposure::TINYINT   AS exposure,
    row_number() OVER () AS user_id
FROM read_csv_auto(getvariable('gz_path'), header = true);
