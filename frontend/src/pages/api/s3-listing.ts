import type { NextApiRequest, NextApiResponse } from "next";
import HTTP_STATUS_CODE from "src/common/constants/HTTP_STATUS_CODE";

interface S3File {
  key: string;
  size: number;
  lastModified: string;
}

type ResponseData = {
  files?: S3File[];
  error?: string;
};

const ALLOWED_BUCKETS = ["gene-expression-assets-public-prod"];
const DEFAULT_REGION = "us-west-2";

function parseS3Response(xml: string): S3File[] {
  const files: S3File[] = [];
  const contentsRegex = /<Contents>([\s\S]*?)<\/Contents>/g;
  let match;

  while ((match = contentsRegex.exec(xml)) !== null) {
    const content = match[1];

    const keyMatch = content.match(/<Key>(.*?)<\/Key>/);
    const sizeMatch = content.match(/<Size>(.*?)<\/Size>/);
    const lastModifiedMatch = content.match(
      /<LastModified>(.*?)<\/LastModified>/
    );

    const key = keyMatch?.[1] || "";
    const size = parseInt(sizeMatch?.[1] || "0", 10);
    const lastModified = lastModifiedMatch?.[1] || "";

    // Skip folder entries
    if (key && !key.endsWith("/")) {
      files.push({ key, size, lastModified });
    }
  }

  return files;
}

export default async function handler(
  request: NextApiRequest,
  response: NextApiResponse<ResponseData>
) {
  if (request.method !== "GET") {
    response
      .status(HTTP_STATUS_CODE.METHOD_NOT_ALLOWED)
      .json({ error: "Method not allowed" });
    return;
  }

  const { bucket, prefix = "", region = DEFAULT_REGION } = request.query;

  if (!bucket || typeof bucket !== "string") {
    response
      .status(HTTP_STATUS_CODE.BAD_REQUEST)
      .json({ error: "Missing bucket parameter" });
    return;
  }

  // Security: Only allow specific buckets
  if (!ALLOWED_BUCKETS.includes(bucket)) {
    response
      .status(HTTP_STATUS_CODE.FORBIDDEN)
      .json({ error: "Bucket not allowed" });
    return;
  }

  const prefixParam =
    typeof prefix === "string" && prefix
      ? `&prefix=${encodeURIComponent(prefix)}`
      : "";

  const regionParam = typeof region === "string" ? region : DEFAULT_REGION;
  const s3Url = `https://${bucket}.s3.${regionParam}.amazonaws.com/?list-type=2${prefixParam}`;

  try {
    const s3Response = await fetch(s3Url);

    if (!s3Response.ok) {
      response.status(HTTP_STATUS_CODE.BAD_GATEWAY).json({
        error: `S3 returned ${s3Response.status}: ${s3Response.statusText}`,
      });
      return;
    }

    const xml = await s3Response.text();
    const files = parseS3Response(xml);

    response.status(HTTP_STATUS_CODE.OK).json({ files });
  } catch (error) {
    response.status(HTTP_STATUS_CODE.INTERNAL_SERVER_ERROR).json({
      error: error instanceof Error ? error.message : "Failed to fetch from S3",
    });
  }
}
