import {
  Popover as PopoverRaw,
  PopoverInteractionKind,
  Position,
  Tag,
} from "@blueprintjs/core";
import { FC } from "react";
import { LeftAlignedDetailsCell } from "../Row/common/style";
import { ContentColumn, ContentWrapper, FieldValues } from "./style";

interface Props {
  values: string[];
}
const CHUNK_SIZE = 25;
const Popover: FC<Props> = ({ values }) => {
  const chunkedValues = Array(Math.ceil(values.length / CHUNK_SIZE))
    .fill("")
    .map((_, index) => index * CHUNK_SIZE)
    .map((begin) => values.slice(begin, begin + CHUNK_SIZE));
  return (
    <LeftAlignedDetailsCell>
      <FieldValues>
        {values[0]}
        <br />
        {values[1]}
      </FieldValues>
      {values.length > 2 && (
        <PopoverRaw
          interactionKind={PopoverInteractionKind.HOVER}
          position={Position.RIGHT}
          boundary="window"
          modifiers={{
            hide: { enabled: false },
            preventOverflow: { enabled: false },
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
          <Tag minimal>+{values.length - 2}</Tag>
        </PopoverRaw>
      )}
    </LeftAlignedDetailsCell>
  );
};

export default Popover;
