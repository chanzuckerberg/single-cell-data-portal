import React from "react";

import { Dialog, InputRadio, DialogTitle } from "@czi-sds/components";
import { useConnect } from "./connect";
import { EmbeddingButtonProps } from "./types";
import CopyButton from "src/components/Collections/components/Dataset/components/DownloadDataset/components/Content/components/DownloadLink/components/CopyButton";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import Highlight from "react-highlight";
import { RadioGroup } from "@mui/material";
import { StyledDialogContent, Label, CodeSnippet, Break } from "./style";
import { StyledButton } from "../../style";

function EmbeddingButton(props: EmbeddingButtonProps) {
  const { project, uniqueMetadata } = props;
  const {
    isOpen,
    language,
    codeSnippet,
    projectTier,
    uri,
    notebookLink,
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
              <InputRadio
                disabled={projectTier !== "maintained"}
                label={projectTier === "maintained" ? "R" : "R (coming soon)"}
                value="r"
              />
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
                    ...uniqueMetadata,
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
                        ...uniqueMetadata,
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
              <a
                href={notebookLink}
                onClick={() =>
                  track(EVENTS.CENSUS_EMBEDDING_NOTEBOOK_CLICKED, {
                    project: project.title,
                    category: projectTier,
                    version: language,
                    ...uniqueMetadata,
                  })
                }
              >
                {language === "python" ? "Jupyter" : "R"} Notebook
              </a>
              !
            </div>
          </div>
        </StyledDialogContent>
      </Dialog>
    </>
  );
}

export default EmbeddingButton;
