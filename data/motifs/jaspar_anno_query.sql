.header on
.mode csv

SELECT
  (m.BASE_ID || '.' || m.VERSION || '.' || m.NAME) AS TfName,
  replace(m.NAME, '::', '+') as GeneName,
  GROUP_CONCAT(DISTINCT CASE WHEN a.TAG = 'family' THEN a.VAL END) AS Family,
  replace(replace(GROUP_CONCAT(DISTINCT p.ACC),',',';'),'-','.') AS UniProt
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
