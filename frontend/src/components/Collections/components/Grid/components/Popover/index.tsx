import {
  Popover as PopoverRaw,
  PopoverInteractionKind,
  Position,
  Tag,
} from "@blueprintjs/core";
import { FC } from "react";
import { PluralizedMetadataLabel } from "src/common/constants/metadata";
import { LeftAlignedDetailsCell } from "../Row/common/style";
import { ContentColumn, ContentWrapper, FieldValues } from "./style";

interface Props {
  label: PluralizedMetadataLabel;
  values: string[];
}
const CHUNK_SIZE = 25;
const Popover: FC<Props> = ({ label, values }) => {
  const chunkedValues = Array(Math.ceil(values.length / CHUNK_SIZE))
    .fill("")
    .map((_, index) => index * CHUNK_SIZE)
    .map((begin) => values.slice(begin, begin + CHUNK_SIZE));
  return (
    <LeftAlignedDetailsCell>
      {values.length <= 2 ? (
        <FieldValues>
          {values[0]}
          <br />
          {values[1]}
        </FieldValues>
      ) : (
        <PopoverRaw
          interactionKind={PopoverInteractionKind.HOVER}
          placement={Position.RIGHT}
          boundary="viewport"
          modifiers={{
            hide: { enabled: false },
            preventOverflow: { enabled: true },
          }}
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
        >
          <Tag minimal>
            {values.length} {label}
          </Tag>
        </PopoverRaw>
      )}
    </LeftAlignedDetailsCell>
  );
};

export default Popover;
