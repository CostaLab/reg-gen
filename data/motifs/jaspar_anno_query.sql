.header on
.mode csv

SELECT
  m.ID,
  (m.BASE_ID || '.' || m.VERSION || '.' || m.NAME) AS tf,
  GROUP_CONCAT(DISTINCT p.ACC) AS UniProt_AC,
  (GROUP_CONCAT(DISTINCT CASE WHEN a.TAG = 'symbol' THEN a.VAL END) || ';' ||
  GROUP_CONCAT(DISTINCT CASE WHEN a.TAG = 'remap_tf_name' THEN a.VAL END)) as name,
  GROUP_CONCAT(DISTINCT CASE WHEN a.TAG = 'family' THEN a.VAL END) AS family
FROM MATRIX m
  JOIN (SELECT *
        FROM MATRIX_ANNOTATION a2
        WHERE NOT exists
        (
            SELECT *
            FROM MATRIX_ANNOTATION a3
            WHERE a2.ID = a3.ID AND a3.TAG = 'tax_group' AND a3.VAL <> 'vertebrates'
        )
        ORDER BY ID) a
    ON m.ID = a.ID
  JOIN MATRIX_PROTEIN p
    ON m.ID = p.ID
WHERE m.COLLECTION = 'CORE'
GROUP BY m.ID, m.COLLECTION, m.BASE_ID, m.VERSION, m.NAME
ORDER BY m.ID;
