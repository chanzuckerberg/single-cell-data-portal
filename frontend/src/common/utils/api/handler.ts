import type { NextApiRequest, NextApiResponse } from "next";
import HTTP_STATUS_CODE from "src/common/constants/HTTP_STATUS_CODE";
import { readJson } from "src/common/utils/api/readJson";

const getHandler = (filePath: string) => {
  const data = readJson(filePath);

  function handler(req: NextApiRequest, res: NextApiResponse) {
    try {
      const entityId = req.query.entityId;
      const response = data[entityId as keyof typeof data];

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
