import {
  Popover,
  PopoverInteractionKind,
  Position,
  Tag,
} from "@blueprintjs/core";
import { PluralizedMetadataLabel } from "src/common/constants/metadata";
import { ContentColumn, ContentWrapper, FieldValues } from "./style";

const CHUNK_SIZE = 25;

interface Props {
  label: PluralizedMetadataLabel;
  values: string[];
}

export default function NTag({ label, values }: Props): JSX.Element {
  const chunkedValues = Array(Math.ceil(values.length / CHUNK_SIZE))
    .fill("")
    .map((_, index) => index * CHUNK_SIZE)
    .map((begin) => values.slice(begin, begin + CHUNK_SIZE));
  return (
    <Popover
      boundary="viewport"
      content={
        <ContentWrapper>
          {chunkedValues.map((chunk, index) => (
            <ContentColumn key={index}>
              {chunk.map((val, idx) => (
                <FieldValues key={val}>
                  {val}
                  {idx !== chunk.length - 1 && <br />}
                </FieldValues>
              ))}
            </ContentColumn>
          ))}
        </ContentWrapper>
      }
      interactionKind={PopoverInteractionKind.HOVER}
      modifiers={{
        hide: { enabled: false },
        preventOverflow: { enabled: true },
      }}
      placement={Position.RIGHT}
    >
      <Tag minimal>
        {values.length} {label}
      </Tag>
    </Popover>
  );
}
