import {
  Popover as PopoverRaw,
  PopoverInteractionKind,
  Position,
  Tag,
} from "@blueprintjs/core";
import React, { FC } from "react";
import { LeftAlignedDetailsCell } from "../common/style";
import { ContentWrapper, FieldValues } from "./style";

interface Props {
  values: string[];
}

const Popover: FC<Props> = ({ values }) => {
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
          position={Position.BOTTOM}
          boundary="window"
          modifiers={{
            hide: { enabled: false },
            preventOverflow: { enabled: false },
          }}
          content={
            <ContentWrapper>
              {values.map((val, idx) => (
                <FieldValues key={val}>
                  {val}
                  {idx !== values.length - 1 && <br />}
                </FieldValues>
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
