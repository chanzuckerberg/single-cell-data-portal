NEW

The initial strate of the state machine is as follows:

- DatasetUploadStatus.WAITING
- DatasetProcessingStatus.INITIALIZED

```mermaid
---
title: Cellxgene Ingestion State Machine
---
stateDiagram-v2
    Validate: Validate
    ValidateAnndata: ValidateAnndata
    note left of ValidateAnndata
        DatasetStatusKey.VALIDATION, DatasetValidationStatus.VALIDATING
        1 download anndata file
        DatasetStatusKey.H5AD, DatasetValidationStatus.VALIDATING
        2 validate anndata file
        DatasetStatusKey.H5AD, DatasetValidationStatus.VALID
        3 upload anndata file
        Fail: DatasetStatusKey.H5AD, DatasetValidationStatus.INVALID
            DatasetStatusKey.VALIDATION, DatasetProcessingStatus.FAILURE
    end note
    ValidateFragment: ValidateFragment
    note left of ValidateFragment
        1 download fragment file and anndata file
        DatasetStatusKey.ATAC_SEQ_FRAGMENT, DatasetValidationStatus.VALIDATING
        2 validate fragment file
        DatasetStatusKey.ATAC_SEQ_FRAGMENT, DatasetValidationStatus.VALID
        DatasetStatusKey.ATAC_SEQ_FRAGMENT, DatasetConversionStatus.CONVERTING
        3 generate index
        DatasetStatusKey.ATAC_SEQ_FRAGMENT, DatasetConversionStatus.CONVERTED
        DatasetStatusKey.ATAC_SEQ_FRAGMENT, DatasetConversionStatus.UPLOADING
        4 upload fragment file and index
        DatasetStatusKey.ATAC_SEQ_FRAGMENT, DatasetUploadStatus.UPLOADED
        Fail: DatasetStatusKey.FRAGMENT, DatasetValidationStatus.INVALID
    end note
    AddLabels: AddLabels
    note left of AddLabels
        1 download original h5ad
        DatasetStatusKey.H5AD, DatasetConversionStatus.CONVERTING
        1 Add labels to anndata file
        DatasetStatusKey.H5AD, DatasetConversionStatus.CONVERTED
        2 extract metadata
        DatasetStatusKey.H5AD, DatasetConversionStatus.UPLOADING
        3 upload anndata file
        DatasetStatusKey.H5AD, DatasetConversionStatus.UPLOADED
        DatasetStatusKey.VALIDATION, DatasetProcessingStatus.SUCCESS
        Fail: DatasetStatusKey.H5AD, DatasetConversionStatus.FAILED
              DatasetStatusKey.H5AD, DatasetProcessingStatus.H5AD
    end note
    CxgSeuratParallel: CxgSeuratParallel
    HandleSuccess: HandleSuccess
    HandleErrors: HandleErrors
    CheckForErrors: CheckForErrors
    ConversionError: ConversionError
    DownloadValidateError: DownloadValidateError
    EndPass: EndPass
    Cxg: Cxg
    CatchCxgFailure: CatchCxgFailure

    [*] --> Validate: manifest
    state Validate {
        [*] --> ValidateAnndata: manifest
        ValidateAnndata --> [*]
        state hasFragment <<choice>>
        [*] --> hasFragment: manifest
        hasFragment --> ValidateFragment: yes
        hasFragment --> [*]: no
        ValidateFragment --> [*]
    }
    Validate --> AddLabels
    AddLabels --> Cxg
    Cxg --> CatchCxgFailure
    Cxg --> HandleSuccess
    HandleSuccess --> CheckForErrors
    CheckForErrors --> DownloadValidateError: $.error
    CheckForErrors --> ConversionError: $[0].error or $[1].error
    CheckForErrors --> EndPass: Default
    Validate --> HandleErrors
    AddLabels --> HandleErrors
    HandleErrors --> CheckForErrors
    ConversionError --> [*]
    DownloadValidateError --> [*]
    EndPass --> [*]
```

Useful States:

- ingestion Initialized
- ingest
- Validating
- Invalid
-
