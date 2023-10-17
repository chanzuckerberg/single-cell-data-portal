import type { NextApiRequest, NextApiResponse } from "next";
import HTTP_STATUS_CODE from "src/common/constants/HTTP_STATUS_CODE";

type ResponseData = {
  "Data Portal"?: string;
  error?: string;
};

export default function handler(
  _: NextApiRequest,
  response: NextApiResponse<ResponseData>
) {
  const commit = process.env.COMMIT_SHA || "";

  if (!commit) {
    response.status(HTTP_STATUS_CODE.INTERNAL_SERVER_ERROR).json({
      error: "Commit SHA is not defined",
    });
  }

  response.status(HTTP_STATUS_CODE.OK).json({
    "Data Portal": commit,
  });
}
