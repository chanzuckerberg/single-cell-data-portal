import React from "react";
import { StyledButton } from "../../../../style";

import { Dialog, InputRadio, DialogTitle } from "@czi-sds/components";
import { useConnect } from "./connect";
import { EmbeddingButtonProps } from "./types";
import CopyButton from "src/components/Collections/components/Dataset/components/DownloadDataset/components/Content/components/DownloadLink/components/CopyButton";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import Highlight from "react-highlight";
import { RadioGroup } from "@mui/material";
import Link from "next/link";
import { StyledDialogContent, Label, CodeSnippet, Break } from "./style";

function EmbeddingButton(props: EmbeddingButtonProps) {
  const { project } = props;
  const {
    isOpen,
    language,
    codeSnippet,
    projectTier,
    uri,
    codeSnippetRef,
    uriTopPosition,
    lineHeight,
    setLanguage,
    handleButtonClick,
  } = useConnect(props);

  if (project.tier === "community") return null;

  return (
    <>
      <StyledButton
        sdsType="primary"
        sdsStyle="square"
        onClick={handleButtonClick}
      >
        Embedding
      </StyledButton>
      <Dialog open={isOpen} onClose={handleButtonClick}>
        <DialogTitle title="Embedding" onClose={handleButtonClick} />
        <StyledDialogContent>
          <div>
            <Label>Language</Label>
            <RadioGroup
              value={language}
              onChange={(_, value: string) => setLanguage(value)}
              row
            >
              <InputRadio label="Python" value="python" />
              <InputRadio disabled label="R (coming soon)" value="r" />
            </RadioGroup>
          </div>
          <div>
            <Label>Quick Start</Label>
            <CodeSnippet
              ref={codeSnippetRef}
              uriTopPosition={uriTopPosition}
              lineHeight={lineHeight}
            >
              <Highlight className={language}>{codeSnippet}</Highlight>
              <CopyButton
                downloadLink={codeSnippet}
                handleAnalytics={() =>
                  track(EVENTS.CENSUS_EMBEDDING_COPIED, {
                    project: project.title,
                    category: projectTier,
                    version: language,
                  })
                }
              />
              {projectTier === "hosted" && (
                <div>
                  <CopyButton
                    downloadLink={uri}
                    handleAnalytics={() =>
                      track(EVENTS.CENSUS_EMBEDDING_COPIED, {
                        project: project.title,
                        category: projectTier,
                        version: "URI",
                      })
                    }
                  />
                </div>
              )}
            </CodeSnippet>
          </div>
          <Break />
          <div>
            <Label>Code Examples</Label>
            <div>
              {
                "If you'd like to see more advanced access patterns, explore this "
              }
              <Link
                href={
                  projectTier === "maintained"
                    ? "https://chanzuckerberg.github.io/cellxgene-census/notebooks/api_demo/census_access_maintained_embeddings.html"
                    : "https://chanzuckerberg.github.io/cellxgene-census/notebooks/api_demo/census_embedding.html"
                }
                onClick={() =>
                  track(EVENTS.CENSUS_EMBEDDING_NOTEBOOK_CLICKED, {
                    project: project.title,
                    category: projectTier,
                    version: language,
                  })
                }
              >
                Jupyter Notebook
              </Link>
              !
            </div>
          </div>
        </StyledDialogContent>
      </Dialog>
    </>
  );
}

export default EmbeddingButton;
