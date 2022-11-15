import logger 

def log_batch_environment():
    batch_environment_variables = [
        "AWS_BATCH_CE_NAME",
        "AWS_BATCH_JOB_ATTEMPT",
        "AWS_BATCH_JOB_ID",
        "STEP_NAME",
        "DROPBOX_URL",
        "ARTIFACT_BUCKET",
        "CELLXGENE_BUCKET",
        "DATASET_ID",
        "DEPLOYMENT_STAGE",
        "MAX_ATTEMPTS",
    ]
    env_vars = dict()
    for var in batch_environment_variables:
        env_vars[var] = os.getenv(var)
    logger.info(f"Batch Job Info: {env_vars}")


def main():
    log_batch_environment()
    dataset_id = os.environ["DATASET_ID"]
    step_name = os.environ["STEP_NAME"]
    return_value = 0
    logger.info(f"Processing dataset {dataset_id}")
    try:
        if step_name == "download-validate":
            from backend.corpora.dataset_processing.process_download_validate import (
                process,
            )

            process(dataset_id, os.environ["DROPBOX_URL"], os.environ["ARTIFACT_BUCKET"])
        elif step_name == "cxg":
            from backend.corpora.dataset_processing.process_cxg import process

            process(
                dataset_id,
                os.environ["ARTIFACT_BUCKET"],
                os.environ["CELLXGENE_BUCKET"],
            )
        elif step_name == "seurat":
            from backend.corpora.dataset_processing.process_seurat import process

            process(dataset_id, os.environ["ARTIFACT_BUCKET"])
        elif step_name == "cxg_remaster":
            try:
                from backend.corpora.dataset_processing.remaster_cxg import process

                process(dataset_id, os.environ["CELLXGENE_BUCKET"], dry_run=False)
            except Exception as e:
                print(f"Dataset {dataset_id} failed to remaster: {e}")
        else:
            logger.error(f"Step function configuration error: Unexpected STEP_NAME '{step_name}'")

    except ProcessingCancelled:
        cancel_dataset(dataset_id)
    except (ValidationFailed, ProcessingFailed, ConversionFailed) as e:
        (status,) = e.args
        update_db(dataset_id, processing_status=status)
        logger.exception("An Error occurred while processing.")
        return_value = 1
    except Exception as e:
        logger.exception(f"An unexpected error occurred while processing the data set: {e}")
        if step_name == "download-validate":
            update_db(
                dataset_id,
                processing_status={"upload_status": UploadStatus.FAILED, "upload_message": str(e)},
            )
        elif step_name == "seurat":
            update_db(dataset_id, processing_status={"rds_status": ConversionStatus.FAILED})
        elif step_name == "cxg":
            update_db(dataset_id, processing_status={"cxg_status": ConversionStatus.FAILED})
        return_value = 1

    return return_value

if __name__ == "__main__":
    rv = main()
    sys.exit(rv)
