import type { NextApiRequest, NextApiResponse } from "next";
import HTTP_STATUS_CODE from "src/common/constants/HTTP_STATUS_CODE";
import { readJson } from "src/common/utils/api/readJson";

const getHandler = (fileOrFolderPath: string) => {
  const isFilePath = fileOrFolderPath.endsWith(".json");
  let data: any;
  if (isFilePath) {
    data = readJson(fileOrFolderPath);
  }

  function handler(req: NextApiRequest, res: NextApiResponse) {
    try {
      const entityId = req.query.entityId;
      let response: any;
      if (!isFilePath && !data) {
        /**
         * (alec): This is a temporary hack to get around the fact that
         * the computational marker genes are stored as individual JSON files where
         * each cell type (entityId) corresponds to a particular file.
         * The entire NextJS API will be removed anyway once the CellGuide pipeline
         * is finished and all CellGuide data fixtures are moved to the cloud.
         */
        const entityIdStr = entityId as string;

        // (alec) This is okay to do on the fly as it only takes up to 5 ms based on empirical testing
        response = readJson(
          `${fileOrFolderPath}/${entityIdStr.replace(":", "_")}.json`
        );
      } else {
        response = data[entityId as keyof typeof data];
      }

      if (response) {
        res.status(HTTP_STATUS_CODE.OK).json(response);
      } else {
        /**
         * (thuang): Setting .status() is not enough. We still need
         * to send the response with .send() or .json().
         * Otherwise, the response will hang and the client will
         * wait forever with a pending request.
         */
        res.status(HTTP_STATUS_CODE.NO_CONTENT).send(null);
      }
    } catch (error) {
      console.error(`Error fetching data: ${error}`);
      res
        .status(HTTP_STATUS_CODE.INTERNAL_SERVER_ERROR)
        .json({ error: "Error fetching data" });
    }
  }

  return handler;
};

export default getHandler;
