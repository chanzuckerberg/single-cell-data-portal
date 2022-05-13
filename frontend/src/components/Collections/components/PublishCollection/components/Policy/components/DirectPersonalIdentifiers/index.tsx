import {
  Popover,
  PopoverInteractionKind,
  PopoverPosition,
} from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import { FC } from "react";
import { GRAY } from "src/components/common/theme";
import {
  BulletWrapper,
  ContentWrapper,
  FirstSentence,
  StyledIcon,
  Text,
} from "./style";

const DirectIdentifiers: FC = ({ children }) => {
  return (
    <Popover
      minimal
      interactionKind={PopoverInteractionKind.HOVER}
      content={<Content />}
      position={PopoverPosition.AUTO}
    >
      {children}
    </Popover>
  );
};

const CATEGORIES = [
  "names",
  "For US: all geographic subdivisions smaller than a state, except for the initial three digits of the ZIP code if the geographic unit formed by combining all ZIP codes with the same three initial digits contains more than 20,000 people",
  "For non-US: any geographic subdivision smaller than a city (e.g. neighborhood, mailing address).",
  "all elements of dates except year, and all ages over 89 or elements indicative of such age",
  "telephone numbers ",
  "fax numbers",
  "email addresses ",
  "government ID numbers (e.g. social security)",
  "medical record numbers",
  "health plan beneficiary numbers",
  "account numbers",
  "certificate or license numbers",
  "vehicle identifiers and license plate numbers",
  "device identifiers and serial numbers",
  "URLs",
  "IP addresses",
  "biometric identifiers",
  "full-face photographs and any comparable images",
  <span key="last">
    or any other unique, identifying characteristic or code, except as permitted
    for re-identification in the{" "}
    <a
      target="_blank"
      rel="noopener"
      href="https://www.hhs.gov/hipaa/for-professionals/privacy/index.html"
    >
      Privacy Rule
    </a>
    .
  </span>,
];

function Content() {
  return (
    <ContentWrapper>
      <div>
        <FirstSentence>
          <Text>
            The below categories are considered “direct personal identifiers”
            pursuant to the{" "}
            <a
              target="_blank"
              rel="noopener"
              href="https://www.hhs.gov/hipaa/for-professionals/privacy/index.html"
            >
              HIPAA Privacy Rule{" "}
            </a>
            as well as other regulations (e.g. GDPR):
          </Text>
        </FirstSentence>
      </div>
      <div>
        {CATEGORIES.map((category, index) => {
          return <Bullet key={index} content={category} />;
        })}
      </div>
    </ContentWrapper>
  );
}

function Bullet({ content }: { content: string | JSX.Element }) {
  return (
    <BulletWrapper>
      <StyledIcon color={GRAY.A} iconSize={4} icon={IconNames.SYMBOL_CIRCLE} />
      <Text>{content}</Text>
    </BulletWrapper>
  );
}

export default DirectIdentifiers;
