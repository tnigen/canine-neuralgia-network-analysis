# Computational methods

## Input processing

The program reads the first worksheet of `drug_target_dog.xlsx` using pandas and openpyxl. It requires the columns `drug` and `dog_target_protein` and refuses an input that already contains `string_protein_id`.

Missing gene-symbol values become empty lookup values. Other symbols are converted to strings and trimmed in a separate representation. The original input values, duplicate relationships, row order, and column order are retained. Unique non-empty lookup symbols are collected in their first-occurrence order.

## STRING queries

Symbols are sent by POST to `https://string-db.org/api/tsv/get_string_ids`. Each unique symbol is submitted in one batch of at most 100. Identifiers are separated by carriage returns. Parameters are `species=9615`, `echo_query=1`, and `caller_identity=dog_target_string_mapper`.

The connection timeout is 10 seconds and the response timeout is 60 seconds. There is a one-second pause between batches. Failed requests are recorded as errors for every gene in the affected batch, and execution continues with subsequent batches. No retry is performed.

## Response parsing and selection

TSV responses must contain `queryIndex`, `stringId`, `ncbiTaxonId`, and `preferredName`. A valid zero-based `queryIndex` links a returned row to its submitted symbol. If that index is invalid, the parser attempts an exact match against the echoed `queryItem`. A row that cannot be linked raises a response error for the batch.

For each symbol, mapping selection proceeds in this order:

1. An empty candidate list is classified as `unmapped`.
2. Candidates are restricted to taxonomy ID `9615`. Absence of canine candidates is an `error`.
3. Canine IDs must start with `9615.ENSCAFP`. Repeated IDs are deduplicated in response order. Absence of a valid ID is an `error`.
4. One unique valid ID is accepted as `mapped`.
5. If several valid IDs remain, one candidate with a preferred name exactly matching the trimmed input symbol is accepted.
6. Remaining ambiguity is classified as `error`, leaving the STRING ID blank.

HTTP 404 classifies every symbol in the batch as `unmapped`. Other request or response errors classify the batch as `error`.

## Outputs and counts

The main output inserts `string_protein_id` immediately after `dog_target_protein`. Empty input symbols receive blank STRING IDs and do not appear in the unique-symbol audit tables.

The QC workbook contains each unique non-empty lookup symbol, its selected STRING ID, and mapping status. The unresolved workbook contains all symbols whose status is not `mapped`, with error or non-match details.

The console mapping rate is `100 × mapped unique symbols / all unique non-empty symbols`, displayed to one decimal place. Its `Unmapped` count includes both unresolved symbols and errors; errors are also reported separately.

## Archived run

The archived workbooks retain their original bytes. The available counts and provenance limits are documented in [RUN_PROVENANCE.md](RUN_PROVENANCE.md). Neither the historical STRING release nor the original API responses are available in this bundle.
