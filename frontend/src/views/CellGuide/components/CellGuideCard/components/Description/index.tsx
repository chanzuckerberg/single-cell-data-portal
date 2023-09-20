import React, { Dispatch, SetStateAction, useEffect, useState } from "react";
import { Icon, TooltipProps } from "@czi-sds/components";
import {
  CellGuideCardDescription,
  ChatGptTooltipSubtext,
  ChatGptTooltipText,
  DescriptionHeader,
  Source,
  SourceLink,
  StyledTooltip,
  Wrapper,
  DescriptionWrapper,
  StyledButton,
  StyledSynonyms,
  FlexContainer,
  StyledOntologyId,
  ValidatedWrapper,
  StyledTag,
  ReferencesWrapper,
  ValidatedInlineWrapper,
} from "./style";
import { Label } from "src/components/Synonyms/style";

import {
  useGptDescription,
  useCellTypeMetadata,
  useValidatedDescription,
} from "src/common/queries/cellGuide";
import Link from "../common/Link";
import { StyledLink } from "../common/Link/style";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { CELL_GUIDE_CORRECTION_SURVEY_LINK } from "src/common/constants/airtableLinks";
import questionMarkIcon from "src/common/images/question-mark-icon.svg";
import { StyledIconImage } from "../common/HelpTooltip/style";
import {
  CELL_GUIDE_CARD_CL_DESCRIPTION,
  CELL_GUIDE_CARD_GPT_DESCRIPTION,
  CELL_GUIDE_CARD_GPT_TOOLTIP_LINK,
  DESCRIPTION_BREAKPOINT_HEIGHT_PX,
  CELL_GUIDE_CARD_SYNONYMS,
} from "src/views/CellGuide/components/CellGuideCard/components/Description/constants";
import { useIsComponentPastBreakpointHeight } from "../common/hooks/useIsComponentPastBreakpoint";

const SLOT_PROPS: TooltipProps["slotProps"] = {
  tooltip: {
    style: {
      maxWidth: 550, // This is needed because SDS bug where width prop doesn't affect dark sdsStyle
      textAlign: "start", // Also dark sdsStyle aligns content to center by default
    },
  },
};

interface DescriptionProps {
  cellTypeName: string;
  cellTypeId: string;
  skinnyMode: boolean;
  setTooltipContent: Dispatch<
    SetStateAction<{
      title: string;
      element: JSX.Element;
    } | null>
  >;
  inSideBar?: boolean;
  synonyms?: string[];
}
export default function Description({
  cellTypeId,
  cellTypeName,
  skinnyMode,
  inSideBar,
  setTooltipContent,
  synonyms,
}: DescriptionProps): JSX.Element {
  const cellTypeIdRaw = cellTypeId.replace(":", "_");
  const [descriptionGpt, setDescriptionGpt] = useState<string>("");
  const [descriptionValidated, setDescriptionValidated] = useState<string>("");
  const [descriptionValidatedReferences, setDescriptionValidatedReferences] =
    useState<string[]>([]);
  const [descriptionCl, setDescriptionCl] = useState<string>("");
  const [descriptionMaxHeight, setDescriptionMaxHeight] = useState<
    number | undefined
  >(DESCRIPTION_BREAKPOINT_HEIGHT_PX);

  const [timerId, setTimerId] = useState<NodeJS.Timer | null>(null); // For chatgpt hover event
  const { isPastBreakpoint, containerRef } = useIsComponentPastBreakpointHeight(
    DESCRIPTION_BREAKPOINT_HEIGHT_PX
  );

  useEffect(() => {
    if (isPastBreakpoint) {
      setDescriptionMaxHeight(DESCRIPTION_BREAKPOINT_HEIGHT_PX);
    } else {
      setDescriptionMaxHeight(undefined);
    }
  }, [isPastBreakpoint]);

  const { data: rawDescriptionGpt } = useGptDescription(cellTypeId);
  const { data: rawDescriptionValidated, isLoading } =
    useValidatedDescription(cellTypeId);
  const { data: cellTypesById } = useCellTypeMetadata();
  const rawDescriptionCl = cellTypesById?.[cellTypeId].clDescription;

  useEffect(() => {
    if (rawDescriptionValidated) {
      setDescriptionValidated(rawDescriptionValidated.description);
      setDescriptionValidatedReferences(rawDescriptionValidated.references);
    } else {
      setDescriptionValidated("");
      setDescriptionValidatedReferences([]);
    }

    if (rawDescriptionGpt) setDescriptionGpt(rawDescriptionGpt);
    else setDescriptionGpt("");

    if (rawDescriptionCl) setDescriptionCl(rawDescriptionCl);
    else setDescriptionCl("");
  }, [rawDescriptionGpt, rawDescriptionCl, rawDescriptionValidated]);

  const copyHandler = () => {
    if (!window) return;

    const selectedText = window.getSelection()?.toString().trim();

    if (selectedText !== "") {
      track(EVENTS.CG_COPY_CELL_TYPE_DESCRIPTION);
    }
  };

  const tooltipContent = (
    <div>
      <ChatGptTooltipText>
        {`This summary on "${cellTypeName}" was generated with ChatGPT, powered by the GPT4 model. Keep in mind that ChatGPT may occasionally present information that is not entirely accurate. For transparency, the prompts used to generate this summary are shared below. CZI is currently offering this as a pilot feature and we may update or change this feature at our discretion.`}
      </ChatGptTooltipText>
      <br />
      <ChatGptTooltipSubtext>
        System role: You are a knowledgeable cell biologist that has
        professional experience writing and curating accurate and informative
        descriptions of cell types.
        <br />
        <br />
        {`User role: I am making a knowledge-base about cell types. Each cell type is a term from the Cell Ontology and will have its own page with a detailed description of that cell type and its function. Please write me a description for "${cellTypeName}". Please return only the description and no other dialogue. The description should include information about the cell type's function. The description should be at least three paragraphs long.`}
      </ChatGptTooltipSubtext>
    </div>
  );

  const submitCorrection = (
    <>
      If you believe a description is inaccurate, please{" "}
      <a
        href={CELL_GUIDE_CORRECTION_SURVEY_LINK}
        target="_blank"
        rel="noreferrer noopener"
        onClick={() => {
          track(EVENTS.SUBMIT_CORRECTION_CLICKED, {
            cell_type_name: cellTypeName,
          });
        }}
      >
        submit a correction
      </a>
      .
    </>
  );
  const disclaimerMessage = (
    <div>
      <em>
        We&apos;re still validating ChatGPT descriptions with our Biocurator
        team.
        {!isPastBreakpoint &&
          " Once a description is validated, we'll add references and a validation icon."}{" "}
        {submitCorrection}
      </em>
    </div>
  );

  const sourceLink = (
    <SourceLink>
      {"Source: ChatGPT "}
      <StyledTooltip
        sdsStyle="dark"
        leaveDelay={0}
        placement="left"
        width="wide"
        arrow
        slotProps={SLOT_PROPS}
        title={!skinnyMode && tooltipContent}
      >
        <StyledLink
          data-testid={CELL_GUIDE_CARD_GPT_TOOLTIP_LINK}
          onMouseOver={() => {
            const id = setTimeout(() => {
              track(EVENTS.CG_CHAT_GPT_HOVER);
            }, 2 * 1000);
            setTimerId(id);
          }}
          onMouseOut={() => {
            if (timerId) {
              clearTimeout(timerId);
              setTimerId(null);
            }
          }}
          onClick={() => {
            skinnyMode &&
              setTooltipContent({
                title: "ChatGPT Descriptions",
                element: tooltipContent,
              });
          }}
          // only setting mobile Tooltip View on touch and not click
          onTouchEnd={() => {
            setTooltipContent({
              title: "ChatGPT Descriptions",
              element: tooltipContent,
            });
          }}
        >
          <StyledIconImage src={questionMarkIcon} />
        </StyledLink>
      </StyledTooltip>
    </SourceLink>
  );
  const readMoreOrLessComponent = (
    <StyledButton
      sdsType="primary"
      sdsStyle="minimal"
      onClick={() => {
        descriptionMaxHeight
          ? setDescriptionMaxHeight(undefined)
          : setDescriptionMaxHeight(DESCRIPTION_BREAKPOINT_HEIGHT_PX);
        descriptionMaxHeight && track(EVENTS.CG_DESCRIPTION_READ_MORE_CLICKED);
      }}
    >
      {descriptionMaxHeight ? "Read More" : "Read Less"}
    </StyledButton>
  );
  const sourceComponent = (
    <Source>
      {isPastBreakpoint ? readMoreOrLessComponent : disclaimerMessage}
      {sourceLink}
    </Source>
  );

  const gptDescriptionComponent = (
    <CellGuideCardDescription
      data-testid={CELL_GUIDE_CARD_GPT_DESCRIPTION}
      onCopy={copyHandler}
      inSideBar={inSideBar}
    >
      <DescriptionWrapper
        inSideBar={inSideBar}
        maxHeight={isPastBreakpoint ? descriptionMaxHeight : undefined}
      >
        {!inSideBar && (
          <DescriptionHeader>Experimental Description</DescriptionHeader>
        )}
        <div ref={containerRef}>
          {inSideBar ? descriptionGpt.split("\n").at(0) : descriptionGpt}
        </div>
      </DescriptionWrapper>
      {!inSideBar && sourceComponent}
      {isPastBreakpoint && !inSideBar && <Source>{disclaimerMessage}</Source>}
    </CellGuideCardDescription>
  );

  const clDescriptionComponent = (
    <>
      <CellGuideCardDescription
        data-testid={CELL_GUIDE_CARD_CL_DESCRIPTION}
        onCopy={copyHandler}
      >
        {descriptionCl}
        <Source>
          <SourceLink>
            {"Source: "}
            <Link
              url={`https://www.ebi.ac.uk/ols4/ontologies/cl/classes/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252F${cellTypeId.replace(
                ":",
                "_"
              )}`}
              label={"Cell Ontology"}
            />
          </SourceLink>
        </Source>
      </CellGuideCardDescription>
      <br />
    </>
  );

  const footerComponent = (
    <FlexContainer>
      {skinnyMode && (
        <StyledOntologyId
          url={`https://www.ebi.ac.uk/ols4/ontologies/cl/classes/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252F${cellTypeIdRaw}`}
          ontologyId={cellTypeId}
        />
      )}
      <StyledSynonyms
        synonyms={synonyms}
        data-testid={CELL_GUIDE_CARD_SYNONYMS}
      />
      {descriptionValidated && (
        <ReferencesWrapper>
          <Label>Citations</Label>
          {descriptionValidatedReferences.map((ref, index) => {
            return (
              <Link
                key={`${ref}-${index}`}
                url={ref}
                label={`[${index + 1}]`}
              />
            );
          })}
        </ReferencesWrapper>
      )}
    </FlexContainer>
  );
  const validatedDescriptionComponent = (
    <CellGuideCardDescription
      data-testid={CELL_GUIDE_CARD_GPT_DESCRIPTION}
      onCopy={copyHandler}
      inSideBar={true} // applies styling from sidebar (no gray background)
    >
      <div>
        <DescriptionWrapper
          inSideBar={inSideBar}
          maxHeight={isPastBreakpoint ? descriptionMaxHeight : undefined}
        >
          <div ref={containerRef}>
            {inSideBar
              ? descriptionValidated.split("\n").at(0)
              : descriptionValidated}
          </div>
        </DescriptionWrapper>
        {isPastBreakpoint && readMoreOrLessComponent}
      </div>
      <ValidatedWrapper>
        <ValidatedInlineWrapper>
          <StyledTag
            sdsType="secondary"
            sdsStyle="square"
            color="success"
            icon={<Icon sdsType="static" sdsIcon="checkCircle" sdsSize="l" />}
            label="Validated"
          />{" "}
          This description has been validated by our Biocurator team.{" "}
          {submitCorrection}
        </ValidatedInlineWrapper>
      </ValidatedWrapper>

      {footerComponent}
    </CellGuideCardDescription>
  );

  const experimentalDescriptionComponent = (
    <>
      {descriptionCl &&
        !inSideBar &&
        !descriptionValidated &&
        clDescriptionComponent}
      {gptDescriptionComponent}
      {footerComponent}
    </>
  );
  const descriptionComponent = descriptionValidated
    ? validatedDescriptionComponent
    : experimentalDescriptionComponent;
  return <Wrapper>{!isLoading && descriptionComponent}</Wrapper>;
}
