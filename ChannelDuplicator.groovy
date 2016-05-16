import groovyx.gpars.dataflow.DataflowQueue

class ChannelDuplicator {
    DataflowQueue orig

    ChannelDuplicator(DataflowQueue orig) {
        this.orig = orig
    }

    DataflowQueue onCopy() {
	    def copied, keep
        (copied, keep) = this.orig.into(2)
        this.orig = keep
        return copied
    }

    static ChannelDuplicator createFrom(Object[] o) {
        return new ChannelDuplicator(Channel.from(o))
    }

    static ChannelDuplicator createFrom(DataflowQueue q) {
        return new ChannelDuplicator(q)
    }
}
