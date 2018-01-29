.header on
.mode csv

SELECT
  trim((m.BASE_ID || '.' || m.VERSION || '.' || m.NAME)) AS TfName,
  trim(replace(m.NAME, '::', '+')) as GeneName,
  trim(replace(GROUP_CONCAT(DISTINCT CASE WHEN a.TAG = 'family' THEN a.VAL END),',',';')) AS Family,
  trim(replace(GROUP_CONCAT(DISTINCT p.ACC),',',';')) AS UniProt,
  trim(GROUP_CONCAT(DISTINCT CASE WHEN a.tag = 'type' THEN a.val END)) AS Source
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
