import { useQuery } from "react-query";
import { API_URL } from "src/configs/configs";
import { API } from "../API";
import { DatasetUploadStatus } from "../entities";
import { apiTemplateToUrl } from "../utils/apiTemplateToUrl";
import { DEFAULT_FETCH_OPTIONS } from "./common";
import { ENTITIES } from "./entities";

export const USE_DATASET_STATUS = {
  entities: [ENTITIES.DATASET_STATUS],
  id: "dataset",
};

async function fetchDatasetStatus(
  _: unknown,
  dataset_uuid: string
): Promise<DatasetUploadStatus> {
  const url = apiTemplateToUrl(API_URL + API.DATASET_STATUS, { dataset_uuid });
  return await (await fetch(url, DEFAULT_FETCH_OPTIONS)).json();
}

export function useDatasetStatus(dataset_uuid: string, shouldFetch: boolean) {
  return useQuery<DatasetUploadStatus>(
    [USE_DATASET_STATUS, dataset_uuid],
    fetchDatasetStatus,
    { enabled: shouldFetch, refetchInterval: 10 * 1000 }
  );
}
